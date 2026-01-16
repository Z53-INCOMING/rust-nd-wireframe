use gif::{Frame, Encoder, Repeat};
use macroquad::miniquad::window::set_window_size;
use macroquad::prelude::*;
use macroquad::rand::rand;
use na::Vector2;
use nalgebra::{self as na, DMatrix, DVector};
use std;
use std::f32::consts::TAU;
use std::fs::File;
use std::fs::read_to_string;

fn rotate_matrix(axis_1: usize, axis_2: usize, angle_in_radians: f32, dimension: usize) -> DMatrix<f32> {
    let mut matrix = DMatrix::identity(dimension, dimension);
    
    matrix[axis_1 + (axis_1 * dimension)] = f32::cos(angle_in_radians);
    matrix[axis_2 + (axis_1 * dimension)] = f32::sin(angle_in_radians);
    
    matrix[axis_1 + (axis_2 * dimension)] = -f32::sin(angle_in_radians);
    matrix[axis_2 + (axis_2 * dimension)] = f32::cos(angle_in_radians);
    
    return matrix;
}

fn load_polytope(path: String, vertices: &mut Vec<DVector<f32>>, edges: &mut Vec<usize>, dimension: &mut usize) {
    if !std::path::Path::new(&path).exists() {
        return
    }

    let contents = std::fs::read_to_string(path).unwrap();
    
    let mut state: u8 = 0;
    
    let mut full_lines_seen = 0;
    
    for line in contents.lines() {
        if line.starts_with("#") {
            continue;
        }
        
        if line.is_empty() {
            if state == 1 { // If done reading rank, start reading vertices
                if full_lines_seen == 1 {
                    state = 2;
                    continue;
                } else {
                    continue;
                }
            } else if state == 2 { // If done reading vertices, start reading edges (faces)
                state = 3;
                continue;
            } else if state == 3 { // If done reading edges (faces), stop
                break;
            }
        }
        
        if line.ends_with("OFF") {
            *dimension = line[.. line.len() - 3].parse().unwrap();
            if *dimension < 4 {
                *dimension = 4;
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
            
            while vertex.len() < 4 {
                vertex.push(0.0);
            }
            
            vertices.push(DVector::from_vec(vertex));
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
        return rotate_matrix(1, axis, (previous_mouse_pos.y - mouse_position().1) * sensitivity, dimension) * rotate_matrix(0, axis, (previous_mouse_pos.x - mouse_position().0) * sensitivity, dimension) * shape_matrix;
    } else {
        return shape_matrix;
    }
}

fn project_vertex(vertex: &DVector<f32>, render_size: f32, screen_size: Vector2<f32>) -> Vector2<f32> {
    let mut screen_vertex = Vector2::new(vertex[0], vertex[1]) / (vertex[2]);
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

fn color_from_off_axis(vertex: &DVector<f32>, w_scale: f32, dimension: usize) -> Color {
    if dimension < 4 {
        return WHITE;
    }
    
    let mut color = Color { r: 0.5, g: 0.5, b: 0.5, a: 1.0 };
    
    let mut distance: f32 = 0.0;
    let mut off_axis: Vec<f32> = vec![];
    for axis in 0..(dimension - 3) {
        distance += vertex[axis + 3] * vertex[axis + 3];
        off_axis.push(vertex[axis + 3]);
    }
    distance = f32::sqrt(distance);
    
    for axis in 0..(dimension - 3) {
        off_axis[axis] /= distance * 2.0;
        off_axis[axis] *= w_scale;
        off_axis[axis] += 0.5;
    }
    
    if dimension > 3 {
        color.r = off_axis[0];
    }
    if dimension > 4 {
        color.g = off_axis[1];
    }
    if dimension > 5 {
        color.b = off_axis[2];
    }
    // if dimension > 6 {
    //     color.a *= off_axis[3];
    // }
    
    color
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
            
            let radius_1 = ((screen_size.y * edge_width) / vertex_1[2]);
            let radius_2 = ((screen_size.y * edge_width) / vertex_2[2]);
            
            let edge_center = vertex_a.lerp(&vertex_b, (s as f32) / ((subdivisions - 1) as f32));
            
            let mut color = if shape_matrix.ncols() < 4 {WHITE} else {color_from_wv(&edge_center, w_scale)};
            // let mut color = color_from_off_axis(&edge_center, w_scale, dimension);
            color.a *= fade_from_depth(edge_center[2], near, far, zoom);
            color.a *= 1.0 - (distance_from_nvolume(&edge_center, 5) * w_scale).clamp(0.0, 1.0);
            
            draw_variable_width_line(project_vertex(&vertex_1, render_size, screen_size), project_vertex(&vertex_2, render_size, screen_size), radius_1 * render_size, radius_2 * render_size, color);
        }
        
    }
    
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
    
    load_polytope("./7 op 5-orthoplex.off".to_string(), &mut vertices, &mut edges, &mut dimension);
    
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
    
    let mut subdivisions = 6;
    
    let mut image_index = -2;
    let frame_count = 120;
    
    set_window_size(1024, 1024);

    loop {
        // Rotate Shape
        if is_mouse_button_pressed(MouseButton::Left) || is_mouse_button_pressed(MouseButton::Middle) {
            (previous_mouse_pos.x, previous_mouse_pos.y) = mouse_position();
        }
        
        if is_mouse_button_down(MouseButton::Left) || is_mouse_button_down(MouseButton::Middle) {
            if is_key_down(KeyCode::Z) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 4, -1.0/216.0);
            } else if is_key_down(KeyCode::X) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 5, -1.0/216.0);
            } else if is_key_down(KeyCode::C) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 6, -1.0/216.0);
            } else if is_key_down(KeyCode::V) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 7, -1.0/216.0);
            } else if is_key_down(KeyCode::LeftControl) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 3, -1.0/216.0);
            } else {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 2, 1.0/216.0);
            }
            
            (previous_mouse_pos.x, previous_mouse_pos.y) = mouse_position();
        }
        
        let scroll = mouse_wheel().1;
        if scroll < 0.0 {
            if is_key_down(KeyCode::LeftControl) {
                render_size *= 12.0/13.0;
            } else if is_key_down(KeyCode::LeftShift) {
                edge_width *= 12.0/13.0;
            } else {
                zoom *= 13.0/12.0;
            }
        } else if scroll > 0.0 {
            if is_key_down(KeyCode::LeftControl) {
                render_size *= 13.0/12.0;
            } else if is_key_down(KeyCode::LeftShift) {
                edge_width *= 13.0/12.0;
            } else {
                zoom *= 12.0/13.0;
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
            shape_matrix = rotate_matrix(0, 2, TAU / (frame_count as f32), shape_matrix.ncols()) * shape_matrix;
            shape_matrix = rotate_matrix(1, 4, TAU / (frame_count as f32), shape_matrix.ncols()) * shape_matrix;
            // shape_matrix = rotate_matrix(1, 3, TAU / (frame_count as f32), shape_matrix.ncols()) * shape_matrix;
            
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
        }
        
        if is_key_pressed(KeyCode::Enter) {
            image_index = -1;
            set_window_size(360, 360);
        }
        
        // println!("{}", shape_position[2]);
        
        next_frame().await
    }
}
