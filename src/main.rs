use macroquad::prelude::*;
use na::{Matrix4, Vector4, Vector2};
use nalgebra as na;
use std;
use std::f32::consts::TAU;

#[macroquad::main("MyGame")]
async fn main() {
    const vertices: [Vector4<f32>; 4] = [
        Vector4::new(1.0, 0.0, 0.0, 0.0),
        Vector4::new(0.0, -1.0, 0.0, 0.0),
        Vector4::new(-1.0, 0.0, 0.0, 0.0),
        Vector4::new(0.0, 1.0, 0.0, 0.0),
    ];

    let edges = [
        0, 1,
        1, 2,
        2, 3,
        3, 0,
    ];
    
    let mut shape_matrix = Matrix4::new_scaling(1.0);
    let mut shape_position = Vector4::new(0.0, 0.0, 0.0, 0.0);
    
    let mut camera_matrix = Matrix4::new_scaling(1.0);
    let mut camera_position = Vector4::new(0.0, 0.0, -3.0, 0.0);
    
    let rotation_matrix = Matrix4::new(
        f32::cos(TAU/256.0), f32::sin(TAU/256.0), 0.0, 0.0,
        -f32::sin(TAU/256.0), f32::cos(TAU/256.0), 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    );

    loop {
        shape_matrix = rotation_matrix * shape_matrix;
        
        clear_background(BLACK);
        
        let inv_camera_matrix = camera_matrix.try_inverse().unwrap();
        
        let mut projected_vertices: [Vector2<f32>; vertices.len()] = [Vector2::default(); vertices.len()];

        let mut index = 0;
        for vertex in vertices {
            // Vertex in world space
            let mut transformed_vertex = (shape_matrix * vertex) + shape_position;
            
            // Vertex in camera space
            transformed_vertex -= camera_position;
            transformed_vertex = inv_camera_matrix * transformed_vertex;
            
            // Vertex in screen space
            let mut screen_vertex = Vector2::new(transformed_vertex.x, transformed_vertex.y) / transformed_vertex.z;
            screen_vertex *= screen_height();
            screen_vertex += Vector2::new(screen_width(), screen_height()) / 2.0;
            
            // Store vertex result
            projected_vertices[index] = screen_vertex;
            index += 1;
        }
        
        for i in (0..edges.len()).step_by(2) {
            draw_line(projected_vertices[edges[i]].x, projected_vertices[edges[i]].y, projected_vertices[edges[i+1]].x, projected_vertices[edges[i+1]].y, 3.0, WHITE);
        }

        next_frame().await
    }
}
