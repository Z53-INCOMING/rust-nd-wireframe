use macroquad::prelude::*;
use na::Vector2;
use nalgebra::dvector;
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
    
    for line in contents.lines() {
        if line.starts_with("#") {
            continue;
        }
        
        if line.is_empty() {
            if state == 1 {
                state = 2;
                continue;
            } else if state == 2 {
                state = 3;
                continue;
            } else if state == 3 {
                break;
            }
        }
        
        if line.ends_with("OFF") {
            *dimension = line[.. line.len() - 3].parse().unwrap();
            state = 1;
            continue;
        }
        
        // Vertices
        if state == 2 {
            let mut vertex: Vec<f32> = vec![];
            
            for coordinate in line.split(" ") {
                vertex.push(coordinate.parse().unwrap());
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

fn color_from_w(w: f32, w_scale: f32) -> Color {
    let positive_w_component = f32::clamp(w * w_scale, 0.0, 1.0);
    let negative_w_component = f32::clamp(-w * w_scale, 0.0, 1.0);
    
    if w > 0.0 {
        return Color::new(1.0, f32::lerp(1.0, 0.5, positive_w_component), f32::lerp(1.0, 0.0, positive_w_component), 1.0 - positive_w_component);
    } else {
        return Color::new(f32::lerp(1.0, 0.0, negative_w_component), f32::lerp(1.0, 0.5, negative_w_component), 1.0, 1.0 - negative_w_component);
    }
}

fn distance_from_volume(vertex: &DVector<f32>) -> f32 {
    let mut distance: f32 = 0.0;
    for axis in 0..(vertex.len() - 3) {
        distance += vertex[axis + 3] * vertex[axis + 3];
    }
    
    f32::sqrt(distance)
}

fn distance_from_4volume(vertex: &DVector<f32>) -> f32 {
    let mut distance: f32 = 0.0;
    for axis in 0..(vertex.len() - 4) {
        distance += vertex[axis + 4] * vertex[axis + 4];
    }
    
    f32::sqrt(distance)
}

fn color_from_off_axis(vertex: &DVector<f32>, w_scale: f32, dimension: usize) -> Color {
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

#[macroquad::main("nD Renderer")]
async fn main() {
    let mut dimension = 0;
    
    let mut vertices: Vec<DVector<f32>> = Vec::new();
    let mut edges: Vec<usize> = Vec::new();
    
    load_polytope("./hexelte.off".to_string(), &mut vertices, &mut edges, &mut dimension);
    
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
    
    let subdivisions = 16;

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
        
        clear_background(BLACK);
        
        let mut local_space_vertices: Vec<DVector<f32>> = Vec::new();
        
        let screen_size = Vector2::new(screen_width(), screen_height());

        for vertex in &vertices {
            // Vertex in world/camera space
            let transformed_vertex = (&shape_matrix * vertex) + &shape_position;
            
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
                
                let mut color = color_from_w(edge_center[3], w_scale);
                // let mut color = color_from_off_axis(&edge_center, w_scale, dimension);
                color.a *= fade_from_depth(edge_center[2], near, far, zoom);
                color.a *= 1.0 - (distance_from_4volume(&edge_center) * w_scale).clamp(0.0, 1.0);
                
                draw_variable_width_line(project_vertex(&vertex_1, render_size, screen_size), project_vertex(&vertex_2, render_size, screen_size), radius_1, radius_2, color);
            }
            
        }
        
        for i in 0..local_space_vertices.len() {
            let coord = project_vertex(&local_space_vertices[i], render_size, screen_size);
            
            let mut color = color_from_w(local_space_vertices[i][3], w_scale);
            // let mut color = color_from_off_axis(&local_space_vertices[i], w_scale, dimension);
            
            color.a *= fade_from_depth(local_space_vertices[i][2], near, far, zoom);
            color.a *= 1.0 - (distance_from_4volume(&local_space_vertices[i]) * w_scale).clamp(0.0, 1.0);
            
            draw_circle(coord.x, coord.y, (screen_size.y * edge_width) / local_space_vertices[i][2], color);
        }

        next_frame().await
    }
}
