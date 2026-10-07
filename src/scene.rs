use std::{fs, io, path::PathBuf};

use nalgebra::{DVector, Vector2};
use serde::Deserialize;
use sfml::graphics::Color;

#[derive(Debug)]
pub struct Scene {
    pub polytopes_folder: PathBuf,
    pub polytope_path: PathBuf,
    pub resolution: u32,
    pub frame_count: i32,
    pub facet_expansion: f32,
    pub facet_expansion_rank: usize,
    pub min_dimension: usize,
    pub dimension: usize,
    pub vertices: Vec<DVector<f32>>,
    pub edges: Vec<usize>,
    pub edge_colors: Vec<Color>,
}

impl Scene {
    /// # Panics
    /// Panics when setup.toml or setup.txt is invalid
    pub fn setup(mut args: impl Iterator<Item = String>) -> Self {
        let given_polytope_path: Option<PathBuf> = args.nth(1).map(Into::into);

        // if you can read setup.toml
        if let Ok(bytes) = fs::read("./setup.toml") {
            let scene_config: SceneConfig =
                toml::from_slice(&bytes).expect("error reading setup.toml");

            // convert SceneConfig -> Scene and return
            let mut scene: Self = scene_config.into();

            if let Some(polytope_path) = given_polytope_path {
                scene.polytope_path = polytope_path;
            }

            return scene;
        }

        let setup_file_contents = match fs::read_to_string("./setup.txt") {
            Ok(v) => v,
            Err(err) => match err.kind() {
                io::ErrorKind::NotFound => panic!("no setup.txt file!!!!"),
                _ => panic!("could not read setup.txt: {err}"),
            },
        };
        let lines: Vec<&str> = setup_file_contents.lines().collect();

        Self {
            polytopes_folder: lines[0].into(),
            polytope_path: given_polytope_path.unwrap_or_else(|| lines[1].into()),
            resolution: lines[2].parse().unwrap(),
            frame_count: lines[3].parse().unwrap(),
            min_dimension: lines[4].parse().unwrap(),
            facet_expansion: lines[5].parse().unwrap(),
            facet_expansion_rank: lines[6].parse::<isize>().unwrap().cast_unsigned(), // converts negative values to super high (integer underflow) ones. necessary for relative to rank values
            dimension: 0,
            vertices: vec![],
            edges: vec![],
            edge_colors: vec![],
        }
    }

    pub fn clear_polytope(&mut self) {
        self.vertices.clear();
        self.edges.clear();
    }

    #[must_use]
    #[allow(dead_code, clippy::cast_precision_loss)]
    pub const fn resolution_vector(&self) -> Vector2<f32> {
        Vector2::new(self.resolution as f32, self.resolution as f32)
    }
}

#[derive(Debug, Deserialize)]
struct SceneConfig {
    polytopes_folder: PathBuf,
    polytope_path: PathBuf,

    min_dimension: usize,

    #[serde(default)]
    animation: AnimationConfig,
    #[serde(default)]
    facet_expansion: FacetExpansionConfig,
}

impl From<SceneConfig> for Scene {
    #[allow(clippy::cast_precision_loss)]
    fn from(value: SceneConfig) -> Self {
        Self {
            polytopes_folder: value.polytopes_folder,
            polytope_path: value.polytope_path,
            resolution: value.animation.resolution,
            frame_count: value.animation.frame_count,
            min_dimension: value.min_dimension,
            facet_expansion: value.facet_expansion.facet_expansion,
            facet_expansion_rank: value.facet_expansion.facet_expansion_rank.cast_unsigned(), // converts negative values to super high (integer underflow) ones. necessary for relative to rank values
            dimension: 0,
            vertices: vec![],
            edges: vec![],
            edge_colors: vec![],
        }
    }
}

#[derive(Debug, Deserialize)]
struct AnimationConfig {
    resolution: u32,
    frame_count: i32,
}

impl Default for AnimationConfig {
    fn default() -> Self {
        Self {
            resolution: 1200,
            frame_count: 80,
        }
    }
}

#[derive(Debug, Deserialize)]
struct FacetExpansionConfig {
    #[serde(default, alias = "amount")]
    facet_expansion: f32,
    #[serde(
        default = "FacetExpansionConfig::default_facet_expansion_rank",
        alias = "rank"
    )]
    facet_expansion_rank: isize,
}

impl FacetExpansionConfig {
    const fn default_facet_expansion_rank() -> isize {
        -1
    }
}

impl Default for FacetExpansionConfig {
    fn default() -> Self {
        Self {
            facet_expansion: 0.0,
            facet_expansion_rank: -1,
        }
    }
}
