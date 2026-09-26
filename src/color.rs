use std::f32::consts::TAU;

use nalgebra::DVector;
use sfml::graphics::{glsl::Vec2, Color};

use crate::math::{inverse_lerp, lerp};

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
#[must_use]
pub fn color_from_hue(hue: f32) -> Color {
    // originally (5. + hue * 6.) / 6. but i simplified it -malki
    let kr = f32::fract(5.0 / 6.0 + hue) * 6.0;
    let kg = f32::fract(3.0 / 6.0 + hue) * 6.0;
    let kb = f32::fract(1.0 / 6.0 + hue) * 6.0;

    let r = 1.0 - f32::min(kr, 4.0 - kr).clamp(0.0, 1.0);
    let g = 1.0 - f32::min(kg, 4.0 - kg).clamp(0.0, 1.0);
    let b = 1.0 - f32::min(kb, 4.0 - kb).clamp(0.0, 1.0);

    Color::rgba((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8, 255)
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
#[must_use]
pub fn color_from_wv(vector: &DVector<f32>, w_scale: f32, edge_color: Color) -> Color {
    if vector.len() < 4 {
        return edge_color;
    }

    let wv_vector = Vec2::new(vector[3], {
        if vector.len() == 4 {
            0.0
        } else {
            vector[4]
        }
    });

    let fade_to_color = color_from_hue((wv_vector.y.atan2(wv_vector.x) / TAU) + 0.5 + (1.0 / 12.0));
    let fade_strength = f32::min(wv_vector.length_sq().sqrt() * w_scale, 1.0);

    Color::rgba(
        lerp(
            f32::from(edge_color.r),
            f32::from(fade_to_color.r),
            (fade_strength * 2.0).min(1.0),
        ) as u8,
        lerp(
            f32::from(edge_color.g),
            f32::from(fade_to_color.g),
            (fade_strength * 2.0).min(1.0),
        ) as u8,
        lerp(
            f32::from(edge_color.b),
            f32::from(fade_to_color.b),
            (fade_strength * 2.0).min(1.0),
        ) as u8,
        ((1.0 - fade_strength) * 255.0) as u8,
    )
}

#[must_use]
pub fn fade_from_depth(z: f32, near: f32, far: f32, zoom: f32) -> f32 {
    1.0 - inverse_lerp(near + zoom, far + zoom, z).clamp(0.0, 1.0)
}
