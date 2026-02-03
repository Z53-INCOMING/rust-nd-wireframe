use macroquad::audio::load_sound;
use macroquad::audio::play_sound_once;
use macroquad::miniquad::window::set_window_size;
use macroquad::prelude::*;
use na::Vector2;
use nalgebra::{self as na, DMatrix, DVector};
use std;
use std::f32::consts::TAU;
use std::vec;

fn rotate_matrix(axis_1: usize, axis_2: usize, angle_in_radians: f32, dimension: usize) -> DMatrix<f32> {
    let mut matrix = DMatrix::identity(dimension, dimension);
    
    matrix[axis_1 + (axis_1 * dimension)] = f32::cos(angle_in_radians);
    matrix[axis_2 + (axis_1 * dimension)] = f32::sin(angle_in_radians);
    
    matrix[axis_1 + (axis_2 * dimension)] = -f32::sin(angle_in_radians);
    matrix[axis_2 + (axis_2 * dimension)] = f32::cos(angle_in_radians);
    
    return matrix;
}

fn get_vertices_from_element(polytope_data: &Vec<Vec<Vec<usize>>>, element_vertices: &mut Vec<usize>, element_edges: &mut Vec<usize>, rank: usize, index: usize) {
    // loop through the facet
    let mut i = 0;
    for sub_element in &polytope_data[rank - 2][index] {
        if rank == 2 { // Faces, add vertices
            element_vertices.push(*sub_element);
            
            // Add the edge, which will be a pair of two integers, each pointing to a vertex ID in the global polytope.
            element_edges.push(*sub_element);
            
            // Get the next vertex ID ID, and append it. Remember that polygons are circular.
            let next_vertex_id_id = (i + 1) % polytope_data[rank - 2][index].len();
            element_edges.push(polytope_data[rank - 2][index][next_vertex_id_id]);
        } else { // Non faces, check sub elements
            let mut sub_vertices: Vec<usize> = vec![];
            let mut sub_edges: Vec<usize> = vec![];
            
            get_vertices_from_element(polytope_data, &mut sub_vertices, &mut sub_edges, rank - 1, *sub_element);
            
            // Merge faces and other elements correctly
            for vertex in sub_vertices.iter() {
                if !element_vertices.contains(vertex) {
                    element_vertices.push(vertex.clone());
                }
            }
            for edge in (0..sub_edges.len()).step_by(2) {
                let vertex_index_a = sub_edges[edge];
                let vertex_index_b = sub_edges[edge + 1];
                
                let mut found_duplicate = false;
                for edge_start_index in (0..element_edges.len()).step_by(2) {
                    let edge_start = element_edges[edge_start_index];
                    let edge_end = element_edges[edge_start_index + 1];
                    
                    if (vertex_index_a == edge_start && vertex_index_b == edge_end) || (vertex_index_a == edge_end && vertex_index_b == edge_start) {
                        found_duplicate = true;
                        break;
                    }
                }
                
                if !found_duplicate {
                    element_edges.push(vertex_index_a);
                    element_edges.push(vertex_index_b);
                }
            }
        }
        
        i += 1;
    }
}

fn load_polytope(path: String, vertices: &mut Vec<DVector<f32>>, edges: &mut Vec<usize>, dimension: &mut usize, min_dimension: usize, facet_expansion: f32) {
    if !std::path::Path::new(&path).exists() {
        return
    }

    let contents = std::fs::read_to_string(path).unwrap();
    
    let mut state: u8 = 0;
    
    let mut rank: u8 = 0;
    
    let mut full_lines_seen = 0;
    
    // vertices
    let mut polytope_vertices: Vec<DVector<f32>> = vec![];
    // rank, element, indices referencing previous rank
    let mut polytope_data: Vec<Vec<Vec<usize>>> = vec![];
    
    for line in contents.lines() {
        if line.starts_with("#") {
            continue;
        }
        
        if line.is_empty() {
            if state == 1 { // If done reading rank, start reading vertices
                if full_lines_seen == 1 {
                    state = 2;
                }
            } else if state == 2 { // If done reading vertices, start reading edges (faces)
                state = 3;
                polytope_data.push(vec![]);
            } else if state > 2 { // If done reading edges (faces), continue or stop depending on facet_expansion
                if facet_expansion == 0.0 {
                    break;
                } else {
                    state += 1;
                    if state > rank {
                        break;
                    } else {
                        polytope_data.push(vec![]);
                    }
                }
            }
            
            continue;
        }
        
        if line.ends_with("OFF") {
            if line == "OFF" {
                *dimension = 3; // some dumbasses think that not having a number for 3D is okay, well, it's NOT,
                // it means I have to take time out of MY day to add an edge case for it every single time I
                // make an OFF importer. AGH.
            } else {
                *dimension = line[.. line.len() - 3].parse().unwrap();
            }
            
            rank = *dimension as u8;
            
            if *dimension < min_dimension {
                *dimension = min_dimension;
            }
            
            state = 1;
            continue;
        }
        
        full_lines_seen += 1;
        
        // Vertices
        if state == 2 {
            let mut vertex: Vec<f32> = vec![];
            
            for coordinate in line.split(" ") {
                if !coordinate.is_empty() {
                    vertex.push(coordinate.parse().unwrap());
                }
            }
            
            while vertex.len() < min_dimension {
                vertex.push(0.0);
            }
            
            polytope_vertices.push(DVector::from_vec(vertex));
        }
        
        // Edges (actually faces)
        if state == 3 {
            // stores the vertex indices of the face
            let mut face: Vec<usize> = vec![];
            
            // go through the line of text to find the indices
            let mut index = 0;
            for number_string in line.split(" ") {
                let number: usize = number_string.parse().unwrap();
                
                // the first one is the size of the face. who needs that? I have .len() and I'm not afraid to use it.
                if index != 0 {
                    face.push(number);
                }
                
                index += 1;
            }
            
            polytope_data[0].push(face.clone());
            
            if facet_expansion == 0.0 {
                // loop through the face to get all the edges
                for index in 0..face.len() {
                    let vertex_index_a = face[index];
                    let vertex_index_b = face[(index + 1) % face.len()];
                    
                    // make sure the edge or its opposite aren't in the edges array
                    let mut found_duplicate = false;
                    for edge_start_index in (0..edges.len()).step_by(2) {
                        let edge_start = edges[edge_start_index];
                        let edge_end = edges[edge_start_index + 1];
                        
                        if (vertex_index_a == edge_start && vertex_index_b == edge_end) || (vertex_index_a == edge_end && vertex_index_b == edge_start) {
                            found_duplicate = true;
                            break;
                        }
                    }
                    
                    // add them
                    if !found_duplicate {
                        edges.push(vertex_index_a);
                        edges.push(vertex_index_b);
                    }
                }
            }
        }
        
        if state > 3 {
            // stores the rank n-1 element indices of the rank n element
            let mut element: Vec<usize> = vec![];
            
            // go through the line of text to find the indices
            let mut index = 0;
            for number_string in line.split(" ") {
                let number: usize = number_string.parse().unwrap();
                
                // the first one is the size of the element. who needs that? I have .len() and I'm not afraid to use it.
                if index != 0 {
                    element.push(number);
                }
                
                index += 1;
            }
            
            polytope_data[(state - 3) as usize].push(element);
        }
    }
    
    // we now have the polytope_data.
    if facet_expansion > 0.0 {
        
        // okay, we need to append vertices of every facet, scaled inward towards the average, to the polytope.
        for facet in 0..polytope_data[polytope_data.len() - 1].len() {
            let mut facet_vertices: Vec<usize> = vec![];
            let mut facet_edges: Vec<usize> = vec![];
            
            get_vertices_from_element(&polytope_data, &mut facet_vertices, &mut facet_edges, (rank - 1) as usize, facet);
            
            // once that is done, loop over all the facet_vertices to determine the center.
            let mut facet_center: DVector<f32> = DVector::zeros(*dimension);
            for vertex in facet_vertices.iter() {
                facet_center += &polytope_vertices[*vertex];
            }
            facet_center /= facet_vertices.len() as f32;
            
            // vertices is the vertices of the mesh
            // polytope_vertices is the vertices of the polytope
            // fucking dumbass (for context, this used to say polytope_vertices.len(), and I struggled to find why it wasn't working)
            let past_vertex_count = vertices.len();
            
            // loop over them again, subtracting each one by the center, multiplying by facet_expansion, and then adding the center
            for vertex in facet_vertices.iter() {
                vertices.push(((&polytope_vertices[*vertex] - &facet_center) * facet_expansion) + &facet_center);
            }
            
            for edge in facet_edges.iter() {
                // These variables are very poorly named, so I will explain
                // past_vertex_count is the number of vertices before a facet, so that's our starting point
                // we have edges as references to vertex IDs in the global polytope
                // but we want the edges as reference to vertex IDs in the facet
                // luckily we can just search the facet_vertices array for these IDs and then it'll work
                
                // edge is the global polytope vertex ID to the first or second half of an edge
                // facet_vertices contains the vertex IDs of the facet in relation to the global polytope
                
                edges.push(past_vertex_count + facet_vertices.iter().position(|x| *x == *edge).unwrap());
            }
        }
        
    } else {
        *vertices = polytope_vertices.clone();
    }
}

fn draw_variable_width_line(start_point: Vector2<f32>, end_point: Vector2<f32>, start_radius: f32, end_radius: f32, color: Color) {
    let edge_direction = (end_point - start_point).normalize();
    let left_of_edge = Vector2::new(edge_direction.y, -edge_direction.x);
    let right_of_edge = Vector2::new(-edge_direction.y, edge_direction.x);
    
    draw_triangle(
        vec2(start_point.x + (left_of_edge.x * start_radius), start_point.y + (left_of_edge.y * start_radius)),
        vec2(start_point.x + (right_of_edge.x * start_radius), start_point.y + (right_of_edge.y * start_radius)),
        vec2(end_point.x + (left_of_edge.x * end_radius), end_point.y + (left_of_edge.y * end_radius)),
        color
    );
    
    draw_triangle(
        vec2(end_point.x + (left_of_edge.x * end_radius), end_point.y + (left_of_edge.y * end_radius)),
        vec2(end_point.x + (right_of_edge.x * end_radius), end_point.y + (right_of_edge.y * end_radius)),
        vec2(start_point.x + (right_of_edge.x * start_radius), start_point.y + (right_of_edge.y * start_radius)),
        color
    );
}

fn mouse_control(previous_mouse_pos: Vector2<f32>, dimension: usize, shape_matrix: DMatrix<f32>, axis: usize, sensitivity: f32) -> DMatrix<f32> {
    if axis < dimension {
        return rotate_matrix(1, axis, (mouse_position().1 - previous_mouse_pos.y) * -sensitivity, dimension) * rotate_matrix(0, axis, (mouse_position().0 - previous_mouse_pos.x) * sensitivity, dimension) * shape_matrix;
    } else {
        return shape_matrix;
    }
}

fn project_vertex(vertex: &DVector<f32>, render_size: f32, screen_size: Vector2<f32>) -> Vector2<f32> {
    let mut screen_vertex = Vector2::new(-vertex[0], vertex[1]) / (vertex[2]);
    screen_vertex *= -screen_height() * render_size;
    screen_vertex += screen_size / 2.0;
    
    screen_vertex
}

fn color_from_hue(hue: f32) -> Color {
    let kr = f32::fract((5.0 + hue * 6.0) / 6.0) * 6.0;
    let kg = f32::fract((3.0 + hue * 6.0) / 6.0) * 6.0;
    let kb = f32::fract((1.0 + hue * 6.0) / 6.0) * 6.0;
    
    let r = 1.0 - f32::max(f32::min(f32::min(kr, 4.0 - kr), 1.0), 0.0);
    let g = 1.0 - f32::max(f32::min(f32::min(kg, 4.0 - kg), 1.0), 0.0);
    let b = 1.0 - f32::max(f32::min(f32::min(kb, 4.0 - kb), 1.0), 0.0);
    
    return Color::new(r, g, b, 1.0);
}

fn color_from_wv(vector: &DVector<f32>, w_scale: f32) -> Color {
    if vector.len() < 4 {
        return WHITE;
    }
    
    if vector.len() == 4 {
        let positive_w_component = f32::clamp(vector[3] * w_scale, 0.0, 1.0);
        let negative_w_component = f32::clamp(-vector[3] * w_scale, 0.0, 1.0);

        if vector[3] > 0.0 {
            return Color::new(1.0, f32::lerp(1.0, 0.5, positive_w_component), f32::lerp(1.0, 0.0, positive_w_component), 1.0 - positive_w_component);
        } else {
            return Color::new(f32::lerp(1.0, 0.0, negative_w_component), f32::lerp(1.0, 0.5, negative_w_component), 1.0, 1.0 - negative_w_component);
        }
    }
    
    let wv_vector = Vec2::new(vector[3], vector[4]);
    
    let fade_to_color = color_from_hue((wv_vector.to_angle() / TAU) + 0.5 + (1.0 / 12.0));
    let fade_strength = f32::min(wv_vector.length() * w_scale, 1.0);
    
    return Color::new(
        f32::lerp(1.0, fade_to_color.r, fade_strength),
        f32::lerp(1.0, fade_to_color.g, fade_strength),
        f32::lerp(1.0, fade_to_color.b, fade_strength),
        1.0 - fade_strength
    );
}

fn distance_from_nvolume(vertex: &DVector<f32>, n: usize) -> f32 {
    if vertex.len() < n {
        return 0.0;
    }
    
    let mut distance: f32 = 0.0;
    for axis in 0..(vertex.len() - n) {
        distance += vertex[axis + n] * vertex[axis + n];
    }
    
    f32::sqrt(distance)
}

fn fade_from_depth(z: f32, near: f32, far: f32, zoom: f32) -> f32 {
    1.0 - clamp(f32::inverse_lerp(near + zoom, far + zoom, z), 0.0, 1.0)
}

fn render(vertices: &Vec<DVector<f32>>, edges: &Vec<usize>, subdivisions: i32, shape_matrix: &DMatrix<f32>, shape_position: &DVector<f32>, edge_width: f32, near: f32, far: f32, zoom: f32, w_scale: f32, render_size: f32) {
    clear_background(BLACK);
    
    let mut local_space_vertices: Vec<DVector<f32>> = Vec::new();
    
    let screen_size = Vector2::new(screen_width(), screen_height());

    for vertex in vertices {
        // Vertex in world/camera space
        let transformed_vertex = (shape_matrix * vertex) + shape_position;
        
        // Store vertex result
        local_space_vertices.push(transformed_vertex);
    }
    
    for i in (0..edges.len()).step_by(2) {
        // A and B are the ends of the edges, 1 and 2 are the ends of the sub edges
        let vertex_a = &local_space_vertices[edges[i]];
        let vertex_b = &local_space_vertices[edges[i + 1]];
        
        for s in 0..subdivisions {
            let vertex_1 = vertex_a.lerp(&vertex_b, (s as f32) / (subdivisions as f32));
            let vertex_2 = vertex_a.lerp(&vertex_b, ((s + 1) as f32) / (subdivisions as f32));
            
            let radius_1 = (screen_size.y * edge_width) / vertex_1[2];
            let radius_2 = (screen_size.y * edge_width) / vertex_2[2];
            
            let edge_center = (&vertex_1 + &vertex_2) / 2.0;
            
            let mut color = if shape_matrix.ncols() < 4 {WHITE} else {color_from_wv(&edge_center, w_scale)};
            // let mut color = color_from_off_axis(&edge_center, w_scale, dimension);
            color.a *= fade_from_depth(edge_center[2], near, far, zoom);
            color.a *= 1.0 - (distance_from_nvolume(&edge_center, 5) * w_scale).clamp(0.0, 1.0);
            
            draw_variable_width_line(project_vertex(&vertex_1, render_size, screen_size), project_vertex(&vertex_2, render_size, screen_size), radius_1 * render_size, radius_2 * render_size, color);
        }
        
    }
    
    // Render vertices
    // for i in 0..local_space_vertices.len() {
    //     let coord = project_vertex(&local_space_vertices[i], render_size, screen_size);
        
    //     let mut color = color_from_wv(&local_space_vertices[i], w_scale);
    //     // let mut color = color_from_off_axis(&local_space_vertices[i], w_scale, dimension);
        
    //     color.a *= fade_from_depth(local_space_vertices[i][2], near, far, zoom);
    //     color.a *= 1.0 - (distance_from_nvolume(&local_space_vertices[i], 5) * w_scale).clamp(0.0, 1.0);
        
    //     draw_circle(coord.x, coord.y, (screen_size.y * edge_width * render_size) / local_space_vertices[i][2], color);
    // }
}

#[macroquad::main("nD Renderer")]
async fn main() {
    let mut dimension = 0;
    
    let mut vertices: Vec<DVector<f32>> = Vec::new();
    let mut edges: Vec<usize> = Vec::new();
    
    load_polytope("./5-cell.off".to_string(), &mut vertices, &mut edges, &mut dimension, 4, 0.0);
    
    let mut shape_matrix = DMatrix::identity(dimension, dimension);
    let mut shape_position = DVector::zeros(dimension);
    shape_position[2] = 2.0;
    
    let mut render_size= 0.5;
    let mut edge_width= 1.0 / 84.0;
    let mut zoom = 2.0;
    
    let mut w_scale: f32 = 0.5;
    let mut near = -1.0;
    let mut far = 0.5;
    
    let mut previous_mouse_pos = Vector2::new(0.0, 0.0);
    
    let mut subdivisions = 1;
    
    let mut image_index = -2;
    let frame_count = 240;
    
    let mut rotations: Vec<usize> = vec![];
    let mut rotations_global_vs_local: Vec<bool> = vec![];
    
    let mut starting_position: Vec<f32> = vec![];
    let mut motion: Vec<f32> = vec![];
    
    let done_sound = load_sound("done.wav").await.unwrap();
    
    set_window_size(1024, 1024);

    loop {
        // Rotate Shape
        if is_mouse_button_pressed(MouseButton::Left) || is_mouse_button_pressed(MouseButton::Middle) {
            (previous_mouse_pos.x, previous_mouse_pos.y) = mouse_position();
        }
        
        if is_mouse_button_down(MouseButton::Left) || is_mouse_button_down(MouseButton::Middle) {
            if is_key_down(KeyCode::LeftControl) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 3, -1.0/216.0);
            } else if is_key_down(KeyCode::Z) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 4, -1.0/216.0);
            } else if is_key_down(KeyCode::X) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 5, -1.0/216.0);
            } else if is_key_down(KeyCode::C) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 6, -1.0/216.0);
            } else if is_key_down(KeyCode::V) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 7, -1.0/216.0);
            } else if is_key_down(KeyCode::B) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 8, -1.0/216.0);
            } else if is_key_down(KeyCode::N) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 9, -1.0/216.0);
            } else {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 2, 1.0/216.0);
            }
            
            (previous_mouse_pos.x, previous_mouse_pos.y) = mouse_position();
        }
        
        let scroll = mouse_wheel().1;
        if scroll < 0.0 {
            if is_key_down(KeyCode::LeftControl) {
                zoom *= 13.0/12.0;
                render_size *= 13.0/12.0;
            } else if is_key_down(KeyCode::LeftShift) {
                edge_width *= 12.0/13.0;
            } else {
                render_size *= 12.0/13.0;
            }
        } else if scroll > 0.0 {
            if is_key_down(KeyCode::LeftControl) {
                zoom *= 12.0/13.0;
                render_size *= 12.0/13.0;
            } else if is_key_down(KeyCode::LeftShift) {
                edge_width *= 13.0/12.0;
            } else {
                render_size *= 13.0/12.0;
            }
        }
        shape_position[2] = zoom;
        
        if is_key_down(KeyCode::Q) {
            near += get_frame_time();
        }
        if is_key_down(KeyCode::A) {
            near -= get_frame_time();
        }
        if is_key_down(KeyCode::W) {
            far += get_frame_time();
        }
        if is_key_down(KeyCode::S) {
            far -= get_frame_time();
        }
        if is_key_down(KeyCode::E) {
            w_scale *= 1.0 - get_frame_time();
        }
        if is_key_down(KeyCode::D) {
            w_scale *= 1.0 + get_frame_time();
        }
        if is_key_pressed(KeyCode::R) {
            subdivisions += 1;
        }
        if is_key_pressed(KeyCode::F) {
            subdivisions -= 1;
        }
        
        render(&vertices, &edges, subdivisions, &shape_matrix, &shape_position, edge_width, near, far, zoom, w_scale, render_size);
        
        if image_index > -1 { // During the loop
            for i in (0..rotations.len()).step_by(2) {
                let rotation_matrix = rotate_matrix(rotations[i], rotations[i + 1], TAU / (frame_count as f32), shape_matrix.ncols());
                if rotations_global_vs_local[i / 2] {
                    shape_matrix = &shape_matrix * rotation_matrix;
                } else {
                    shape_matrix = rotation_matrix * &shape_matrix;
                }
            }
            for i in 0..motion.len() {
                if i != 2 {
                    shape_position[i] += motion[i] / (frame_count as f32);
                }
            }
            
            get_screen_data().export_png(&format!("./images/{:03}.png", image_index));
            
            image_index += 1;
        }
        if image_index == -1 {
            image_index = 0;
        }
        
        if image_index == frame_count { // End
            set_default_camera();
            image_index = -2;
            set_window_size(1024, 1024);
            shape_position[0] = 0.0;
            
            play_sound_once(&done_sound);
        }
        
        if is_key_pressed(KeyCode::Enter) { // Start
            image_index = -1;
            
            if !std::path::Path::new("./src/rotations.txt").exists() {
                panic!("no rotations.txt file!!!!");
            }
            if !std::path::Path::new("./src/motion.txt").exists() {
                panic!("no motion.txt file!!!!");
            }

            let rotation_file_contents = std::fs::read_to_string("./src/rotations.txt").unwrap();
            
            rotations.clear();
            rotations_global_vs_local.clear();
            
            let mut rotation_file_values: Vec<usize> = vec![];
            
            for line in rotation_file_contents.lines() {
                // go through the line of text to find the numbers
                for number_string in line.split(" ") {
                    let number: usize = number_string.parse().unwrap();
                    
                    rotation_file_values.push(number);
                }
            }
            
            for i in (0..rotation_file_values.len()).step_by(3) {
                rotations.push(rotation_file_values[i]);
                rotations.push(rotation_file_values[i + 1]);
                
                rotations_global_vs_local.push(rotation_file_values[i + 2] == 1);
            }
            
            let motion_file_contents = std::fs::read_to_string("./src/motion.txt").unwrap();
            
            starting_position.clear();
            motion.clear();
            
            let mut index = 0;
            for line in motion_file_contents.lines() {
                // go through the line of text to find the numbers
                for number_string in line.split(" ") {
                    let number: f32 = number_string.parse().unwrap();
                    
                    if index == 0 {
                        starting_position.push(number);
                    } else if index == 1 {
                        motion.push(number);
                    }
                }
                
                index += 1;
            }
            
            for i in 0..starting_position.len() {
                if i != 2 {
                    shape_position[i] = starting_position[i];
                }
            }
            
            println!("{:?}", starting_position);
            
            // set_window_size(256, 256);
        }
        
        next_frame().await
    }
}
