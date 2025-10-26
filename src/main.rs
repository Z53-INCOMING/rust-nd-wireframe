use macroquad::prelude::*;
use macroquad::rand::rand;
use na::{Matrix4, Vector4, Vector2};
use nalgebra::{self as na, Cholesky, Reflection4};
use std;
use std::f32::consts::TAU;

fn xy_matrix(radians: f32) -> Matrix4<f32> {
    let matrix = Matrix4::new(
        f32::cos(radians), f32::sin(radians), 0.0, 0.0,
        -f32::sin(radians), f32::cos(radians), 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    );
    
    return matrix;
}

fn xz_matrix(radians: f32) -> Matrix4<f32> {
    let matrix = Matrix4::new(
        f32::cos(radians), 0.0, f32::sin(radians), 0.0,
        0.0, 1.0, 0.0, 0.0,
        -f32::sin(radians), 0.0, f32::cos(radians), 0.0,
        0.0, 0.0, 0.0, 1.0
    );
    
    return matrix;
}

fn yz_matrix(radians: f32) -> Matrix4<f32> {
    let matrix = Matrix4::new(
        1.0, 0.0, 0.0, 0.0,
        0.0, f32::cos(radians), f32::sin(radians), 0.0,
        0.0, -f32::sin(radians), f32::cos(radians), 0.0,
        0.0, 0.0, 0.0, 1.0
    );
    
    return matrix;
}

fn yw_matrix(radians: f32) -> Matrix4<f32> {
    let matrix = Matrix4::new(
        1.0, 0.0, 0.0, 0.0,
        0.0, f32::cos(radians), 0.0, f32::sin(radians),
        0.0, 0.0, 1.0, 0.0,
        0.0, -f32::sin(radians), 0.0, f32::cos(radians),
    );
    
    return matrix;
}

fn xw_matrix(radians: f32) -> Matrix4<f32> {
    let matrix = Matrix4::new(
        f32::cos(radians), 0.0, 0.0, f32::sin(radians),
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        -f32::sin(radians), 0.0, 0.0, f32::cos(radians),
    );
    
    return matrix;
}

fn generate_polytope_from_coxeter_diagram(vertices: &mut Vec<Vector4<f32>>, edges: &mut Vec<usize>, coxeter_edges: &Vec<u8>, coxeter_rings: &Vec<bool>) {
    let mut angles: Vec<f32> = Vec::new();
    
    for weight in coxeter_edges {
        // every angle is an nth of a half turn
        angles.push(TAU / ((*weight * 2) as f32));
    }
    
    let mut pre_matrix: Matrix4<f32> = Matrix4::identity();
    
    let rank = coxeter_edges.len() + 1;
    
    for i in 0..rank {
        for j in 0..rank {
            if i != j { // preserve diagonal
                if i8::abs((i as i8) - (j as i8)) == 1 { // if the nodes are right next to each other, there's an edge connecting them
                    pre_matrix[i * 4 + j] = f32::cos(angles[usize::min(i, j)]);
                } // otherwise, the weight is 2, and cosine of 90 degrees is 0, and it's already 0
            }
        }
    }
    
    let mirror_matrix: Matrix4<f32> = pre_matrix.cholesky().expect("invalid coxeter diagram").unpack().transpose();
    
    let mut point: Vector4<f32> = Vector4::new(0.0, 0.0, 0.0, 0.0);
    
    // thanks buzz for the wedge product
    let fundamental_simplex = Matrix4::from_fn(|r, c| {
            (-1f32).powi(r as i32) * mirror_matrix.remove_row(r).remove_column(c).determinant()
        }
    );
    
    let mut ringed_node_count = 0;
    for i in 0..rank {
        if coxeter_rings[i] {
            point += fundamental_simplex.column(i);
            ringed_node_count += 1;
        }
    }
    point /= ringed_node_count as f32;
    point = point.normalize();
    
    for i in 0..256 {
        vertices.push(point);
        
        let mirror_index= (rand() % rank as u32) as usize;
        
        point += -point.dot(&mirror_matrix.column(mirror_index)) * mirror_matrix.column(mirror_index) * 2.0;
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

#[macroquad::main("4D Renderer")]
async fn main() {
    let mut vertices: Vec<Vector4<f32>> = Vec::new();
    let mut edges: Vec<usize> = Vec::new();
    
    generate_polytope_from_coxeter_diagram(&mut vertices, &mut edges, &vec![3, 3, 3], &vec![true, true, true, false]);
    
    let mut shape_matrix = Matrix4::new_scaling(1.0);
    let mut shape_position = Vector4::new(0.0, 0.0, 0.0, 0.0);
    
    let mut camera_matrix = Matrix4::new_scaling(1.0);
    let mut camera_position = Vector4::new(0.0, 0.0, -2.0, 0.0);
    
    let mut render_size= 0.5;
    let mut edge_width= 1.0 / 84.0;
    
    let mut previous_mouse_pos = Vector2::new(0.0, 0.0);

    loop {
        // Rotate Shape
        if is_mouse_button_pressed(MouseButton::Left) || is_mouse_button_pressed(MouseButton::Middle) {
            (previous_mouse_pos.x, previous_mouse_pos.y) = mouse_position();
        }
        
        if is_mouse_button_down(MouseButton::Left) || is_mouse_button_down(MouseButton::Middle) {
            if is_key_down(KeyCode::LeftControl) {
                shape_matrix = xw_matrix((previous_mouse_pos.x - mouse_position().0) / -216.0) * shape_matrix;
                shape_matrix = yw_matrix((previous_mouse_pos.y - mouse_position().1) / -216.0) * shape_matrix;
            } else {
                shape_matrix = xz_matrix((previous_mouse_pos.x - mouse_position().0) / -216.0) * shape_matrix;
                shape_matrix = yz_matrix((previous_mouse_pos.y - mouse_position().1) / -216.0) * shape_matrix;
            }
            
            (previous_mouse_pos.x, previous_mouse_pos.y) = mouse_position();
        }
        
        let scroll = mouse_wheel().1;
        if scroll < 0.0 {
            if is_key_down(KeyCode::LeftControl) {
                render_size *= 12.0/13.0;
            } else {
                camera_position.z *= 13.0/12.0;
            }
        } else if scroll > 0.0 {
            if is_key_down(KeyCode::LeftControl) {
                render_size *= 13.0/12.0;
            } else {
                camera_position.z *= 12.0/13.0;
            }
        }
        
        
        clear_background(BLACK);
        
        let inv_camera_matrix = camera_matrix.try_inverse().unwrap();
        
        let mut projected_vertices: Vec<Vector2<f32>> = Vec::new();
        let mut depth_of_vertices: Vec<f32> = Vec::new();
        
        let screen_size = Vector2::new(screen_width(), screen_height());

        let mut index = 0;
        for vertex in &vertices {
            // Vertex in world space
            let mut transformed_vertex = (shape_matrix * vertex) + shape_position;
            
            // Vertex in camera space
            transformed_vertex -= camera_position;
            transformed_vertex = inv_camera_matrix * transformed_vertex;
            
            // Vertex in screen space
            let mut screen_vertex = Vector2::new(transformed_vertex.x, transformed_vertex.y) / transformed_vertex.z;
            screen_vertex *= -screen_height() * render_size;
            screen_vertex += screen_size / 2.0;
            
            // Store vertex result
            projected_vertices.push(screen_vertex);
            depth_of_vertices.push(transformed_vertex.z);
            index += 1;
        }
        
        for i in (0..edges.len()).step_by(2) {
            let vertex_a = Vector2::new(projected_vertices[edges[i]].x, projected_vertices[edges[i]].y);
            let vertex_b = Vector2::new(projected_vertices[edges[i + 1]].x, projected_vertices[edges[i + 1]].y);
            
            let radius_a = ((screen_size.y * edge_width) / depth_of_vertices[edges[i]]);
            let radius_b = ((screen_size.y * edge_width) / depth_of_vertices[edges[i + 1]]);
            
            draw_variable_width_line(vertex_a, vertex_b, radius_a, radius_b, WHITE);
        }
        
        for i in (0..projected_vertices.len()) {
            draw_circle(projected_vertices[i].x, projected_vertices[i].y, ((screen_size.y * edge_width) / depth_of_vertices[i]), WHITE);
        }

        next_frame().await
    }
}
