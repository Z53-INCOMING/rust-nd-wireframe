use nalgebra::{DMatrix, DVector};
use sfml::cpp::FBox;
use sfml::graphics::glsl::Vec2;
use sfml::graphics::{Color, PrimitiveType, RenderStates, RenderTarget, Vertex};

use crate::color::{color_from_wv, fade_from_depth};
use crate::math::{distance_from_nvolume, normalize, project_vertex};
use crate::scene::Scene;

pub fn draw_triangle_color(
    target: &mut FBox<impl RenderTarget>,
    v1: Vec2,
    v2: Vec2,
    v3: Vec2,
    color1: Color,
    color2: Color,
    color3: Color,
) {
    let vertices = [
        Vertex::with_pos_color(Vec2::new(v1.x, v1.y), color1),
        Vertex::with_pos_color(Vec2::new(v2.x, v2.y), color2),
        Vertex::with_pos_color(Vec2::new(v3.x, v3.y), color3),
    ];

    target.draw_primitives(&vertices, PrimitiveType::TRIANGLES, &RenderStates::DEFAULT);
}

pub fn draw_variable_width_line(
    target: &mut FBox<impl RenderTarget>,
    start_point: Vec2,
    end_point: Vec2,
    start_radius: f32,
    end_radius: f32,
    start_color: Color,
    end_color: Color,
) {
    if start_color.a > 0 || end_color.a > 0 {
        let edge_direction = normalize(end_point - start_point);
        let left_of_edge = Vec2::new(edge_direction.y, -edge_direction.x);
        let right_of_edge = Vec2::new(-edge_direction.y, edge_direction.x);

        draw_triangle_color(
            target,
            start_point + (left_of_edge * start_radius),
            start_point + (right_of_edge * start_radius),
            end_point + (left_of_edge * end_radius),
            start_color,
            start_color,
            end_color,
        );

        draw_triangle_color(
            target,
            end_point + (left_of_edge * end_radius),
            end_point + (right_of_edge * end_radius),
            start_point + (right_of_edge * start_radius),
            end_color,
            end_color,
            start_color,
        );
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FadePlanes {
    /// Fade start
    pub near: f32,
    /// Fade end
    pub far: f32,
    /// Extra-dimensional fade (how far you can see into what's perpendicular to XYZ)
    pub w_scale: f32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeSettings {
    /// How wide the drawn edges are
    pub edge_width: f32,
    /// How much subdivisions the edges have
    pub subdivisions: u8,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CameraPerspective {
    pub zoom: f32,
    /// Multiplier of scale, helps with perspective
    pub render_size: f32,
}

#[allow(clippy::cast_precision_loss)]
pub fn render(
    target: &mut FBox<impl RenderTarget>,
    scene: &Scene,
    shape_matrix: &DMatrix<f32>,
    shape_position: &DVector<f32>,
    edge_settings: EdgeSettings,
    fade_planes: FadePlanes,
    camera: CameraPerspective,
) {
    let EdgeSettings {
        edge_width,
        subdivisions,
    } = edge_settings;
    let FadePlanes { w_scale, .. } = fade_planes;
    let CameraPerspective { render_size, .. } = camera;
    let screen_size: Vec2 = target.size().as_other();

    target.clear(Color::BLACK);

    let mut local_space_vertices: Vec<DVector<f32>> = Vec::new();

    for vertex in &scene.vertices {
        // Vertex in world/camera space
        let transformed_vertex = (shape_matrix * vertex) + shape_position;

        // Store vertex result
        local_space_vertices.push(transformed_vertex);
    }

    for i in (0..scene.edges.len()).step_by(2) {
        // A and B are the ends of the edges, 1 and 2 are the ends of the sub edges
        let vertex_a = &local_space_vertices[scene.edges[i]];
        let vertex_b = &local_space_vertices[scene.edges[i + 1]];

        for s in 0..subdivisions {
            let vertex_1 = vertex_a.lerp(vertex_b, f32::from(s) / f32::from(subdivisions));
            let vertex_2 = vertex_a.lerp(vertex_b, (f32::from(s) + 1.0) / f32::from(subdivisions));

            let radius_1 = (screen_size.y * edge_width) / vertex_1[2];
            let radius_2 = (screen_size.y * edge_width) / vertex_2[2];

            let mut color_1 = color_from_wv(&vertex_1, w_scale, scene.edge_colors[i / 2]);
            color_1.a = get_alpha(color_1, &vertex_1, fade_planes, camera);

            let mut color_2 = color_from_wv(&vertex_2, w_scale, scene.edge_colors[i / 2]);
            color_2.a = get_alpha(color_2, &vertex_2, fade_planes, camera);

            draw_variable_width_line(
                target,
                project_vertex(&vertex_1, render_size, screen_size),
                project_vertex(&vertex_2, render_size, screen_size),
                radius_1 * render_size,
                radius_2 * render_size,
                color_1,
                color_2,
            );
        }
    }
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
#[inline]
fn get_alpha(
    color: Color,
    vertex: &DVector<f32>,
    fade_planes: FadePlanes,
    camera: CameraPerspective,
) -> u8 {
    let mut a = f32::from(color.a) / 255.0;
    a *= fade_from_depth(vertex[2], fade_planes.near, fade_planes.far, camera.zoom);
    a *= 1.0 - (distance_from_nvolume(vertex, 5) * fade_planes.w_scale).clamp(0.0, 1.0);

    (a.clamp(0.0, 1.0) * 255.0) as u8
}
