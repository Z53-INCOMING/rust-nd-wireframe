use macroquad::prelude::*;
use na::Vector2;
use nalgebra::{self as na, DMatrix, DVector};
use std;
use std::f32::consts::TAU;
use std::fs::File;

fn rotate_matrix(axis_1: usize, axis_2: usize, angle_in_radians: f32, dimension: usize) -> DMatrix<f32> {
    let mut matrix = DMatrix::identity(dimension, dimension);
    
    matrix[axis_1 + (axis_1 * dimension)] = f32::cos(angle_in_radians);
    matrix[axis_2 + (axis_1 * dimension)] = f32::sin(angle_in_radians);
    
    matrix[axis_1 + (axis_2 * dimension)] = -f32::sin(angle_in_radians);
    matrix[axis_2 + (axis_2 * dimension)] = f32::cos(angle_in_radians);
    
    return matrix;
}

fn load_polytope(path: String, vertices: &mut Vec<DVector<f32>>, edges: &mut Vec<usize>) {
    
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

fn mouse_control(previous_mouse_pos: Vector2<f32>, dimension: usize, shape_matrix: DMatrix<f32>, axis: usize) -> DMatrix<f32> {
    if axis < dimension {
        return rotate_matrix(1, axis, (previous_mouse_pos.y - mouse_position().1) / 216.0, dimension) * rotate_matrix(0, axis, (previous_mouse_pos.x - mouse_position().0) / 216.0, dimension) * shape_matrix;
    } else {
        return shape_matrix;
    }
}

#[macroquad::main("nD Renderer")]
async fn main() {
    let mut dimension = 5;
    
    let mut vertices: Vec<DVector<f32>> = Vec::new();
    let mut edges: Vec<usize> = Vec::new();
    
    vertices = [
        DVector::from_column_slice(&[-1.0, -1.0, -1.0, -1.0, 0.0]),
        DVector::from_column_slice(&[1.0, -1.0, -1.0, -1.0, 0.0]),
        DVector::from_column_slice(&[-1.0, 1.0, -1.0, -1.0, 0.0]),
        DVector::from_column_slice(&[1.0, 1.0, -1.0, -1.0, 0.0]),
        DVector::from_column_slice(&[-1.0, -1.0, 1.0, -1.0, 0.0]),
        DVector::from_column_slice(&[1.0, -1.0, 1.0, -1.0, 0.0]),
        DVector::from_column_slice(&[-1.0, 1.0, 1.0, -1.0, 0.0]),
        DVector::from_column_slice(&[1.0, 1.0, 1.0, -1.0, 0.0]),
        DVector::from_column_slice(&[-1.0, -1.0, -1.0, 1.0, 0.0]),
        DVector::from_column_slice(&[1.0, -1.0, -1.0, 1.0, 0.0]),
        DVector::from_column_slice(&[-1.0, 1.0, -1.0, 1.0, 0.0]),
        DVector::from_column_slice(&[1.0, 1.0, -1.0, 1.0, 0.0]),
        DVector::from_column_slice(&[-1.0, -1.0, 1.0, 1.0, 0.0]),
        DVector::from_column_slice(&[1.0, -1.0, 1.0, 1.0, 0.0]),
        DVector::from_column_slice(&[-1.0, 1.0, 1.0, 1.0, 0.0]),
        DVector::from_column_slice(&[1.0, 1.0, 1.0, 1.0, 0.0]),
        DVector::from_column_slice(&[-1.0, -1.0, -1.0, -1.0, 1.0])
    ].to_vec();

    edges = [
        0b0000, 0b0001,
        0b0000, 0b0010,
        0b0001, 0b0011,
        0b0010, 0b0011,
        0b0100, 0b0101,
        0b0100, 0b0110,
        0b0101, 0b0111,
        0b0110, 0b0111,
        0b0000, 0b0100,
        0b0001, 0b0101,
        0b0010, 0b0110,
        0b0011, 0b0111,
        0b1000, 0b1001,
        0b1000, 0b1010,
        0b1001, 0b1011,
        0b1010, 0b1011,
        0b1100, 0b1101,
        0b1100, 0b1110,
        0b1101, 0b1111,
        0b1110, 0b1111,
        0b1000, 0b1100,
        0b1001, 0b1101,
        0b1010, 0b1110,
        0b1011, 0b1111,
        0b0000, 0b1000,
        0b0001, 0b1001,
        0b0010, 0b1010,
        0b0011, 0b1011,
        0b0100, 0b1100,
        0b0101, 0b1101,
        0b0110, 0b1110,
        0b0111, 0b1111,
        0b0000, 0b10000
    ].to_vec();
    
    let mut shape_matrix = DMatrix::identity(dimension, dimension);
    let mut shape_position = DVector::zeros(dimension);
    shape_position[2] = 2.0;
    
    let mut render_size= 0.5;
    let mut edge_width= 1.0 / 84.0;
    let mut zoom = 1.0;
    
    let mut previous_mouse_pos = Vector2::new(0.0, 0.0);
    
    let subdivisions = 4;

    loop {
        // Rotate Shape
        if is_mouse_button_pressed(MouseButton::Left) || is_mouse_button_pressed(MouseButton::Middle) {
            (previous_mouse_pos.x, previous_mouse_pos.y) = mouse_position();
        }
        
        if is_mouse_button_down(MouseButton::Left) || is_mouse_button_down(MouseButton::Middle) {
            if is_key_down(KeyCode::Z) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 4);
            } else if is_key_down(KeyCode::X) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 5);
            } else if is_key_down(KeyCode::C) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 6);
            } else if is_key_down(KeyCode::V) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 7);
            } else if is_key_down(KeyCode::LeftControl) {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 3);
            } else {
                shape_matrix = mouse_control(previous_mouse_pos, dimension, shape_matrix, 2);
            }
            
            (previous_mouse_pos.x, previous_mouse_pos.y) = mouse_position();
        }
        
        let scroll = mouse_wheel().1;
        if scroll < 0.0 {
            if is_key_down(KeyCode::LeftControl) {
                render_size *= 12.0/13.0;
            } else {
                shape_position[2] *= 13.0/12.0;
            }
        } else if scroll > 0.0 {
            if is_key_down(KeyCode::LeftControl) {
                render_size *= 13.0/12.0;
            } else {
                shape_position[2] *= 12.0/13.0;
            }
        }
        
        
        clear_background(BLACK);
        
        let mut projected_vertices: Vec<Vector2<f32>> = Vec::new();
        let mut local_space_vertices: Vec<DVector<f32>> = Vec::new();
        
        let screen_size = Vector2::new(screen_width(), screen_height());

        let mut index = 0;
        for vertex in &vertices {
            // Vertex in world/camera space
            let mut transformed_vertex = (&shape_matrix * vertex) + &shape_position;
            
            // Vertex in screen space
            let mut screen_vertex = Vector2::new(transformed_vertex[0], transformed_vertex[1]) / (transformed_vertex[2]);
            screen_vertex *= -screen_height() * render_size;
            screen_vertex += screen_size / 2.0;
            
            // Store vertex result
            projected_vertices.push(screen_vertex);
            local_space_vertices.push(transformed_vertex);
            index += 1;
        }
        
        for i in (0..edges.len()).step_by(2) {
            let vertex_a = Vector2::new(projected_vertices[edges[i]].x, projected_vertices[edges[i]].y);
            let vertex_b = Vector2::new(projected_vertices[edges[i + 1]].x, projected_vertices[edges[i + 1]].y);
            
            let radius_a = ((screen_size.y * edge_width) / local_space_vertices[edges[i]][2]);
            let radius_b = ((screen_size.y * edge_width) / local_space_vertices[edges[i + 1]][2]);
            
            draw_variable_width_line(vertex_a, vertex_b, radius_a, radius_b, WHITE);
        }
        
        for i in (0..projected_vertices.len()) {
            draw_circle(projected_vertices[i].x, projected_vertices[i].y, ((screen_size.y * edge_width) / local_space_vertices[i][2]), WHITE);
        }

        next_frame().await
    }
}
