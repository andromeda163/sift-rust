# sift-rust

Numerically accurate, from-scratch Rust port of VLFeat's SIFT core.

## Current Coverage (numerical layer)

All of the following are implemented and tested in `src/numerical/`:

| Component                | File                        | Tests                                      | C cross-validation      |
|--------------------------|-----------------------------|--------------------------------------------|-------------------------|
| Fast `exp(-x)` approx    | `fast_expn.rs`              | Unit + `test_against_c.sh`                 | Yes (exact match)       |
| Gaussian smoothing       | `smooth.rs`                 | Unit (`vl_sift_smooth`)                    | No (planned)            |
| 2× linear upsampling     | `upsample.rs`               | Unit (`copy_and_upsample_rows`)            | No (planned)            |
| Power-of-2 downsampling  | `downsample.rs`             | Unit (`copy_and_downsample`)               | No (planned)            |

- All public APIs are re-exported from the crate root for compatibility.
- Run `cargo test` for the full suite.
- Run `./test_against_c.sh` for the strict numerical match against VLFeat C (currently only covers `fast_expn`).

## Roadmap (next layers)

1. Scale-space / octave construction (using the above primitives)
2. Difference-of-Gaussians (DoG)
3. Keypoint localization & refinement
4. Orientation assignment
5. Descriptor extraction
6. Full `VlSift` detector + descriptor API + optional FFI

The goal is a complete, numerically verified, dependency-free Rust SIFT that matches VLFeat's output on the same inputs.

## Notes

- This crate is the Rust-native reimplementation effort; the original C code lives in `../vlfeat/`.
- Stage-by-stage binary dump comparison harness exists in the test module for future incremental verification.
