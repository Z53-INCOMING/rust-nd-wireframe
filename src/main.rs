use macroquad::prelude::*;
use na::{Matrix4, Vector4, Vector2};
use nalgebra as na;
use std;
use std::f32::consts::TAU;

const render_size: f32 = 0.5;
const edge_width: f32 = 1.0 / 84.0;

const vertices: [Vector4<f32>; 8] = [
    Vector4::new(1.0, 0.0, 0.0, 0.0),
    Vector4::new(-1.0, 0.0, 0.0, 0.0),
    Vector4::new(0.0, 1.0, 0.0, 0.0),
    Vector4::new(0.0, -1.0, 0.0, 0.0),
    Vector4::new(0.0, 0.0, 1.0, 0.0),
    Vector4::new(0.0, 0.0, -1.0, 0.0),
    Vector4::new(0.0, 0.0, 0.0, 1.0),
    Vector4::new(0.0, 0.0, 0.0, -1.0),
];

const edges: [usize; 48] = [
    0, 2,
    0, 3,
    0, 4,
    0, 5,
    0, 6,
    0, 7,
    1, 2,
    1, 3,
    1, 4,
    1, 5,
    1, 6,
    1, 7,
    2, 4,
    2, 5,
    2, 6,
    2, 7,
    3, 4,
    3, 5,
    3, 6,
    3, 7,
    4, 6,
    4, 7,
    5, 6,
    5, 7,
];

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

fn yw_matrix(radians: f32) -> Matrix4<f32> {
    let matrix = Matrix4::new(
        1.0, 0.0, 0.0, 0.0,
        0.0, f32::cos(radians), 0.0, f32::sin(radians),
        0.0, 0.0, 1.0, 0.0,
        0.0, -f32::sin(radians), 0.0, f32::cos(radians),
    );
    
    return matrix;
}

#[macroquad::main("4D Renderer")]
async fn main() {
    let mut shape_matrix = Matrix4::new_scaling(1.0);
    let mut shape_position = Vector4::new(0.0, 0.0, 0.0, 0.0);
    
    let mut camera_matrix = Matrix4::new_scaling(1.0);
    let mut camera_position = Vector4::new(0.0, 0.0, -2.0, 0.0);

    loop {
        shape_matrix = xy_matrix(TAU * get_frame_time() * 0.125) * shape_matrix;
        shape_matrix = xz_matrix(TAU * get_frame_time() * 0.6165165326523 * 0.25) * shape_matrix;
        shape_matrix = yw_matrix(TAU * get_frame_time() * 0.1440937294359 * 0.25) * shape_matrix;
        
        clear_background(BLACK);
        
        let inv_camera_matrix = camera_matrix.try_inverse().unwrap();
        
        let mut projected_vertices: [Vector2<f32>; vertices.len()] = [Vector2::default(); vertices.len()];
        let mut depth_of_vertices: [f32; vertices.len()] = [0.0; vertices.len()];
        
        let screen_size = Vector2::new(screen_width(), screen_height());

        let mut index = 0;
        for vertex in vertices {
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
            projected_vertices[index] = screen_vertex;
            depth_of_vertices[index] = transformed_vertex.z;
            index += 1;
        }
        
        for i in (0..edges.len()).step_by(2) {
            let vertex_a = Vector2::new(projected_vertices[edges[i]].x, projected_vertices[edges[i]].y);
            let vertex_b = Vector2::new(projected_vertices[edges[i + 1]].x, projected_vertices[edges[i + 1]].y);
            
            let edge_direction = (vertex_b - vertex_a).normalize();
            let left_of_edge = Vector2::new(edge_direction.y, -edge_direction.x);
            let right_of_edge = Vector2::new(-edge_direction.y, edge_direction.x);
            
            let a_radius = ((screen_size.y * edge_width) / depth_of_vertices[edges[i]]);
            let b_radius = ((screen_size.y * edge_width) / depth_of_vertices[edges[i + 1]]);
            
            draw_triangle(
                vec2(vertex_a.x + (left_of_edge.x * a_radius), vertex_a.y + (left_of_edge.y * a_radius)),
                vec2(vertex_a.x + (right_of_edge.x * a_radius), vertex_a.y + (right_of_edge.y * a_radius)),
                vec2(vertex_b.x + (left_of_edge.x * b_radius), vertex_b.y + (left_of_edge.y * b_radius)),
                WHITE
            );
            
            draw_triangle(
                vec2(vertex_b.x + (left_of_edge.x * b_radius), vertex_b.y + (left_of_edge.y * b_radius)),
                vec2(vertex_b.x + (right_of_edge.x * b_radius), vertex_b.y + (right_of_edge.y * b_radius)),
                vec2(vertex_a.x + (right_of_edge.x * a_radius), vertex_a.y + (right_of_edge.y * a_radius)),
                WHITE
            );
        }
        
        for i in (0..projected_vertices.len()) {
            draw_circle(projected_vertices[i].x, projected_vertices[i].y, ((screen_size.y * edge_width) / depth_of_vertices[i]), WHITE);
        }

        next_frame().await
    }
}
