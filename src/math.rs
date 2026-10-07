use nalgebra::{DMatrix, DVector};
use sfml::graphics::glsl::Vec2;

#[inline]
#[must_use]
pub fn lerp(a: f32, b: f32, t: f32) -> f32 {
    t.mul_add(b - a, a)
}

#[inline]
#[must_use]
pub fn inverse_lerp(a: f32, b: f32, v: f32) -> f32 {
    (v - a) / (b - a)
}

#[inline]
#[must_use]
pub fn normalize(p: Vec2) -> Vec2 {
    p / p.length_sq().sqrt()
}

#[must_use]
pub fn rotate_matrix(
    axis_1: usize,
    axis_2: usize,
    angle_in_radians: f32,
    dimension: usize,
) -> DMatrix<f32> {
    let mut matrix = DMatrix::identity(dimension, dimension);

    matrix[axis_1 + (axis_1 * dimension)] = f32::cos(angle_in_radians);
    matrix[axis_2 + (axis_1 * dimension)] = f32::sin(angle_in_radians);

    matrix[axis_1 + (axis_2 * dimension)] = -f32::sin(angle_in_radians);
    matrix[axis_2 + (axis_2 * dimension)] = f32::cos(angle_in_radians);

    matrix
}

#[must_use]
pub fn project_vertex(vertex: &DVector<f32>, render_size: f32, screen_size: Vec2) -> Vec2 {
    let mut screen_vertex = Vec2::new(-vertex[0], vertex[1]) / (vertex[2]);
    screen_vertex *= -screen_size.y * render_size;
    screen_vertex += screen_size / 2.0;

    screen_vertex
}

#[must_use]
pub fn distance_from_nvolume(vertex: &DVector<f32>, n: usize) -> f32 {
    if vertex.len() < n {
        return 0.0;
    }

    let distance: f32 = (0..(vertex.len() - n))
        .map(|axis| vertex[axis + n] * vertex[axis + n])
        .sum();

    f32::sqrt(distance)
}
