use macroquad::prelude::*;
use na::{Matrix4, Vector4};
use nalgebra as na;

#[macroquad::main("MyGame")]
async fn main() {
    let vertices = [
        Vector4::new(1.0, 0.0, 0.0, 0.0),
        Vector4::new(0.0, -1.0, 0.0, 0.0),
        Vector4::new(-1.0, 0.0, 0.0, 0.0),
        Vector4::new(0.0, 1.0, 0.0, 0.0),
    ];

    loop {
        clear_background(BLACK);

        draw_line(64.0, 32.0, 64.0, 96.0, 3.0, WHITE);
        draw_line(128.0, 32.0, 128.0, 96.0, 3.0, WHITE);
        draw_line(160.0, 32.0, 160.0, 96.0, 3.0, WHITE);
        draw_line(64.0, 128.0, 64.0, 128.0 + 64.0, 3.0, WHITE);
        draw_line(96.0, 128.0, 96.0, 128.0 + 64.0, 3.0, WHITE);
        draw_line(128.0, 128.0, 128.0, 128.0 + 64.0, 3.0, WHITE);
        draw_line(140.0, 128.0 + 64.0, 128.0 + 64.0, 128.0 + 64.0, 3.0, WHITE);

        next_frame().await
    }
}
