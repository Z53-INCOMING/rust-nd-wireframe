#![allow(clippy::cast_precision_loss)]

use nalgebra::{DMatrix, DVector, VecStorage};
use sfml::audio::{Sound, SoundBuffer};
use sfml::graphics::glsl::Vec2;
use sfml::graphics::{RenderTarget, RenderTexture, RenderWindow};
use sfml::system::Vector2i;
use sfml::window::mouse::{Button, Wheel};
use sfml::window::{ContextSettings, Event, Key, Style};
use std::collections::HashSet;
use std::env;
use std::f32::consts::TAU;

use crate::loader::load_polytope;
use crate::math::rotate_matrix;
use crate::render::{render, CameraPerspective, EdgeSettings, FadePlanes};
use crate::scene::Scene;

mod color;
mod loader;
mod math;
mod render;
mod scene;

const DONE_SOUND_BYTES: &[u8] = include_bytes!(".././done.wav");
const FPS: u32 = 60;
const FRAME_TIME: f32 = 1.0 / FPS as f32;

fn main() {
    // Create folders if they do not exist already
    std::fs::create_dir_all("./images").unwrap();
    std::fs::create_dir_all("./polytopes").unwrap();

    let mut scene = Scene::setup(env::args());

    load_polytope(&mut scene, false);

    let mut shape_matrix = DMatrix::identity(scene.dimension, scene.dimension);
    let mut shape_position = DVector::zeros(scene.dimension);
    shape_position[2] = 2.0;

    // despite shape_matrix being defined the exact same way, only this variable needs to specify its type. ???
    let mut rotational_offset: nalgebra::Matrix<
        f32,
        nalgebra::Dyn,
        nalgebra::Dyn,
        VecStorage<f32, nalgebra::Dyn, nalgebra::Dyn>,
    > = DMatrix::identity(scene.dimension, scene.dimension);

    let mut fade_planes = FadePlanes {
        near: -1.0,
        far: 0.5,
        w_scale: 0.5,
    };

    let mut edge_settings = EdgeSettings {
        edge_width: 1.0 / 84.0,
        subdivisions: 1,
    };

    let mut camera = CameraPerspective {
        render_size: 0.5,
        zoom: 2.0,
    };

    let mut previous_mouse_pos = Vector2i::new(0, 0);
    let mut mouse_lock: bool = false;
    // 0 = left, 1 = middle, 2 = right
    let mut mouse_buttons: [bool; 3] = [false, false, false];
    let mut mouse_scroll_delta: f32;

    let mut down_keys: HashSet<Key> = HashSet::new();

    let facet_expansion_key_speed = f32::exp2(0.25); // 2 ^ 1/4

    let mut image_index = -2;

    let mut rotations: Vec<usize> = vec![];
    let mut rotations_global_vs_local: Vec<bool> = vec![];
    let mut rotation_amounts: Vec<f32> = vec![];

    let mut starting_position: Vec<f32> = vec![];
    let mut motion: Vec<f32> = vec![];

    let done_sound = SoundBuffer::from_memory(DONE_SOUND_BYTES).expect("invalid done sound");
    let mut sound;

    let mut virtual_image = RenderTexture::new(scene.resolution, scene.resolution).unwrap();
    let mut window = RenderWindow::new(
        (1024, 1024),
        "Rust ND Renderer",
        Style::CLOSE,
        &ContextSettings {
            antialiasing_level: 4,
            ..Default::default()
        },
    )
    .unwrap();
    window.set_framerate_limit(FPS);

    'mainloop: loop {
        mouse_scroll_delta = 0.0;
        let mut enter_pressed = false;
        let mut escape_pressed = false;

        while let Some(ev) = window.poll_event() {
            match ev {
                Event::Closed => break 'mainloop,
                Event::MouseButtonPressed { button, .. } => {
                    let button_index = match button {
                        Button::Left => 0,
                        Button::Middle => 1,
                        Button::Right => 2,
                        _ => continue,
                    };

                    mouse_buttons[button_index] = true;
                }
                Event::MouseButtonReleased { button, .. } => {
                    let button_index = match button {
                        Button::Left => 0,
                        Button::Middle => 1,
                        Button::Right => 2,
                        _ => continue,
                    };

                    mouse_buttons[button_index] = false;
                }
                Event::MouseWheelScrolled {
                    wheel: Wheel::VerticalWheel,
                    delta,
                    ..
                } => mouse_scroll_delta = delta,
                Event::KeyPressed { code, .. } => {
                    down_keys.insert(code);
                    match code {
                        Key::R => {
                            edge_settings.subdivisions =
                                edge_settings.subdivisions.saturating_add(1);
                        }
                        Key::F => {
                            edge_settings.subdivisions =
                                edge_settings.subdivisions.saturating_sub(1);
                        }

                        Key::T => {
                            // increases facet_expansion
                            scene.clear_polytope();
                            scene.facet_expansion =
                                1.0 - (1.0 - scene.facet_expansion) / facet_expansion_key_speed;
                            load_polytope(&mut scene, false);
                        }
                        Key::G => {
                            // decreases facet_expansion
                            scene.clear_polytope();
                            scene.facet_expansion = (1.0 - scene.facet_expansion)
                                .mul_add(-facet_expansion_key_speed, 1.0);
                            if scene.facet_expansion < 1.0 - 1.0 / facet_expansion_key_speed {
                                scene.facet_expansion = 0.0;
                            }
                            load_polytope(&mut scene, false);
                        }
                        Key::Y => {
                            scene.clear_polytope();
                            scene.facet_expansion_rank += 1;
                            if scene.facet_expansion_rank > scene.dimension - 1 {
                                scene.facet_expansion_rank = scene.dimension - 1;
                            }
                            load_polytope(&mut scene, false);
                        }
                        Key::H => {
                            scene.clear_polytope();
                            scene.facet_expansion_rank -= 1;
                            if scene.facet_expansion_rank < 2 {
                                scene.facet_expansion_rank = 2;
                            }
                            load_polytope(&mut scene, false);
                        }

                        Key::Num0 => {
                            scene = Scene::setup(env::args());
                            load_polytope(&mut scene, false);
                            if scene.dimension != shape_position.nrows() {
                                shape_matrix = DMatrix::identity(scene.dimension, scene.dimension);
                                rotational_offset =
                                    DMatrix::identity(scene.dimension, scene.dimension);
                                shape_position = DVector::zeros(scene.dimension);
                            }
                        }
                        Key::K => {
                            scene.clear_polytope();
                            load_polytope(&mut scene, true);
                        }
                        Key::Num1 => {
                            rotational_offset = shape_matrix.clone() * rotational_offset;
                            shape_matrix = DMatrix::identity(scene.dimension, scene.dimension);
                        }
                        Key::J => {
                            mouse_lock = !mouse_lock;
                            window.set_mouse_cursor_visible(!mouse_lock);
                            window.set_mouse_cursor_grabbed(mouse_lock);
                        }
                        Key::Escape => escape_pressed = true,
                        Key::Enter => enter_pressed = true,
                        _ => (),
                    }
                }
                Event::KeyReleased { code, .. } => {
                    down_keys.remove(&code);
                }
                _ => (),
            }
        }

        if mouse_lock || mouse_buttons[0] || mouse_buttons[1] {
            if down_keys.contains(&Key::LControl) {
                shape_matrix = mouse_control(
                    &window,
                    previous_mouse_pos,
                    scene.dimension,
                    shape_matrix,
                    3,
                    -1.0 / 216.0,
                );
            } else if down_keys.contains(&Key::Z) {
                shape_matrix = mouse_control(
                    &window,
                    previous_mouse_pos,
                    scene.dimension,
                    shape_matrix,
                    4,
                    -1.0 / 216.0,
                );
            } else if down_keys.contains(&Key::X) {
                shape_matrix = mouse_control(
                    &window,
                    previous_mouse_pos,
                    scene.dimension,
                    shape_matrix,
                    5,
                    -1.0 / 216.0,
                );
            } else if down_keys.contains(&Key::C) {
                shape_matrix = mouse_control(
                    &window,
                    previous_mouse_pos,
                    scene.dimension,
                    shape_matrix,
                    6,
                    -1.0 / 216.0,
                );
            } else if down_keys.contains(&Key::V) {
                shape_matrix = mouse_control(
                    &window,
                    previous_mouse_pos,
                    scene.dimension,
                    shape_matrix,
                    7,
                    -1.0 / 216.0,
                );
            } else if down_keys.contains(&Key::B) {
                shape_matrix = mouse_control(
                    &window,
                    previous_mouse_pos,
                    scene.dimension,
                    shape_matrix,
                    8,
                    -1.0 / 216.0,
                );
            } else if down_keys.contains(&Key::N) {
                shape_matrix = mouse_control(
                    &window,
                    previous_mouse_pos,
                    scene.dimension,
                    shape_matrix,
                    9,
                    -1.0 / 216.0,
                );
            } else {
                shape_matrix = mouse_control(
                    &window,
                    previous_mouse_pos,
                    scene.dimension,
                    shape_matrix,
                    2,
                    1.0 / 216.0,
                );
            }
        }

        if mouse_buttons[2] {
            let mouse_delta = Vec2::new(
                window.mouse_position().x as f32 - (window.size().x as f32 / 2.0),
                window.mouse_position().y as f32 - (window.size().y as f32 / 2.0),
            ) - (Vec2::new(
                previous_mouse_pos.x as f32 - (window.size().x as f32 / 2.0),
                previous_mouse_pos.y as f32 - (window.size().y as f32 / 2.0),
            ));
            let angle_diff = mouse_delta.y.atan2(mouse_delta.x);

            shape_matrix = rotate_matrix(0, 1, angle_diff, scene.dimension) * &shape_matrix;
        }

        previous_mouse_pos.x = window.mouse_position().x;
        previous_mouse_pos.y = window.mouse_position().y;

        if mouse_scroll_delta < 0.0 {
            if down_keys.contains(&Key::LControl) {
                camera.zoom *= 13.0 / 12.0;
                camera.render_size *= 13.0 / 12.0;
            } else if down_keys.contains(&Key::LShift) {
                edge_settings.edge_width *= 12.0 / 13.0;
            } else {
                camera.render_size *= 12.0 / 13.0;
            }
        } else if mouse_scroll_delta > 0.0 {
            if down_keys.contains(&Key::LControl) {
                camera.zoom *= 12.0 / 13.0;
                camera.render_size *= 12.0 / 13.0;
            } else if down_keys.contains(&Key::LShift) {
                edge_settings.edge_width *= 13.0 / 12.0;
            } else {
                camera.render_size *= 13.0 / 12.0;
            }
        }
        shape_position[2] = camera.zoom;

        if down_keys.contains(&Key::Q) {
            fade_planes.near += FRAME_TIME;
        }
        if down_keys.contains(&Key::A) {
            fade_planes.near -= FRAME_TIME;
        }
        if down_keys.contains(&Key::W) {
            fade_planes.far += FRAME_TIME;
        }
        if down_keys.contains(&Key::S) {
            fade_planes.far -= FRAME_TIME;
        }
        if down_keys.contains(&Key::E) {
            fade_planes.w_scale *= 1.0 - FRAME_TIME;
        }
        if down_keys.contains(&Key::D) {
            fade_planes.w_scale *= 1.0 + FRAME_TIME;
        }

        if image_index > -1 {
            // set camera to render target

            // render the scene to virtual_image
            render(
                &mut virtual_image,
                &scene,
                &(&shape_matrix * &rotational_offset),
                &shape_position,
                edge_settings,
                fade_planes,
                camera,
            );
        }

        // render the scene to the screen
        render(
            &mut window,
            &scene,
            &(&shape_matrix * &rotational_offset),
            &shape_position,
            edge_settings,
            fade_planes,
            camera,
        );

        if image_index > -1 {
            // During the loop
            for i in (0..rotations.len()).step_by(2) {
                let rotation_matrix = rotate_matrix(
                    rotations[i],
                    rotations[i + 1],
                    rotation_amounts[i / 2] / (scene.frame_count as f32),
                    shape_matrix.ncols(),
                );
                if rotations_global_vs_local[i / 2] {
                    shape_matrix = &shape_matrix * rotation_matrix;
                } else {
                    shape_matrix = rotation_matrix * &shape_matrix;
                }
            }
            for i in 0..motion.len() {
                if i != 2 {
                    shape_position[i] += motion[i] / (scene.frame_count as f32);
                }
            }

            virtual_image.display();

            virtual_image
                .texture()
                .copy_to_image()
                .expect("could not copy texture to image")
                .save_to_file(&format!("./images/{image_index:03}.png"))
                .expect("could not save image");

            image_index += 1;
        }
        if image_index == -1 {
            image_index = 0;
        }

        if image_index == scene.frame_count {
            // End
            image_index = -2;
            for i in 0..scene.dimension {
                shape_position[i] = 0.0;
            }

            sound = Sound::new();
            sound.set_buffer(&done_sound);
            sound.play();
        }

        if escape_pressed {
            for i in 0..scene.dimension {
                shape_position[i] = 0.0;
            }
            for i in (0..rotations.len()).step_by(2) {
                let rotation_matrix = rotate_matrix(
                    rotations[i],
                    rotations[i + 1],
                    (rotation_amounts[i / 2] / (scene.frame_count as f32)) * (-image_index) as f32,
                    shape_matrix.ncols(),
                );
                if rotations_global_vs_local[i / 2] {
                    shape_matrix = &shape_matrix * rotation_matrix;
                } else {
                    shape_matrix = rotation_matrix * &shape_matrix;
                }
            }
            image_index = -2;
        }

        if enter_pressed {
            image_index = -1;

            start_animation(
                &mut rotations,
                &mut rotations_global_vs_local,
                &mut rotation_amounts,
                &mut starting_position,
                &mut shape_position,
                &mut motion,
            );
        }

        window.display();
    }
}

fn start_animation(
    rotations: &mut Vec<usize>,
    rotations_global_vs_local: &mut Vec<bool>,
    rotation_amounts: &mut Vec<f32>,
    starting_position: &mut Vec<f32>,
    shape_position: &mut DVector<f32>,
    motion: &mut Vec<f32>,
) {
    // Start
    let rotation_file_contents =
        std::fs::read_to_string("./rotations.txt").expect("could not read rotations.txt!");

    rotations.clear();
    rotations_global_vs_local.clear();
    rotation_amounts.clear();

    let mut rotation_file_values: Vec<usize> = vec![];

    for line in rotation_file_contents.lines() {
        let mut value_count = 0;

        // go through the line of text to find the numbers
        for number_string in line.split(' ') {
            let number: usize = number_string.parse().unwrap();

            rotation_file_values.push(number);
            value_count += 1;
        }

        if value_count == 2 {
            rotation_file_values.push(0);
        }
        if value_count == 3 {
            rotation_file_values.push(1);
        }
    }

    for i in (0..rotation_file_values.len()).step_by(4) {
        rotations.push(rotation_file_values[i]);
        rotations.push(rotation_file_values[i + 1]);

        rotations_global_vs_local.push(rotation_file_values[i + 2] == 1);

        rotation_amounts.push(TAU / (rotation_file_values[i + 3] as f32));
    }

    let motion_file_contents =
        std::fs::read_to_string("./motion.txt").expect("could not read motions.txt!");

    starting_position.clear();
    motion.clear();

    for (index, line) in motion_file_contents.lines().enumerate() {
        // go through the line of text to find the numbers
        for number_string in line.split(' ') {
            let number: f32 = number_string.parse().unwrap();

            if index == 0 {
                starting_position.push(number);
            } else if index == 1 {
                motion.push(number);
            }
        }
    }

    for i in 0..starting_position.len() {
        if i != 2 {
            shape_position[i] = starting_position[i];
        }
    }
}

fn mouse_control(
    window: &RenderWindow,
    previous_mouse_pos: Vector2i,
    dimension: usize,
    shape_matrix: DMatrix<f32>,
    axis: usize,
    sensitivity: f32,
) -> DMatrix<f32> {
    if axis < dimension {
        rotate_matrix(
            1,
            axis,
            (window.mouse_position().y - previous_mouse_pos.y) as f32 * -sensitivity,
            dimension,
        ) * rotate_matrix(
            0,
            axis,
            (window.mouse_position().x - previous_mouse_pos.x) as f32 * sensitivity,
            dimension,
        ) * shape_matrix
    } else {
        shape_matrix
    }
}
