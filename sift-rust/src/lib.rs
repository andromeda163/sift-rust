//! Numerically accurate Rust port of VLFeat SIFT core numerical functions from sift.c.
//! These are the building blocks for the full SIFT detector/descriptor port.

pub mod numerical;

// Re-export the public API at the crate root for compatibility with existing
// callers (example, test_against_c.sh, etc.)
pub use numerical::{copy_and_downsample, copy_and_upsample_rows, fast_expn, fast_expn_init};
pub use numerical::{sift_smooth, vl_sift_smooth};

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Read;
    use approx::assert_abs_diff_eq;
    use float_cmp::approx_eq;

    /// Load a .bin file dumped by the modified VLFeat C code (f32 array).
    fn load_bin(path: &str) -> Option<Vec<f32>> {
        let mut f = File::open(path).ok()?;
        let mut buf = Vec::new();
        f.read_to_end(&mut buf).ok()?;
        // assume nel = file_size / 4
        let nel = buf.len() / 4;
        let mut v = vec![0f32; nel];
        // safety: reinterpret
        unsafe {
            std::ptr::copy_nonoverlapping(buf.as_ptr(), v.as_mut_ptr() as *mut u8, buf.len());
        }
        Some(v)
    }

    /// Compare two slices with both absolute and relative tolerance.
    /// Uses approx + float-cmp for robust stage-by-stage numerical testing.
    fn assert_approx_slices(a: &[f32], b: &[f32], abs_tol: f32, rel_tol: f32) {
        assert_eq!(a.len(), b.len(), "length mismatch in stage dump");
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            assert_abs_diff_eq!(*x, *y, epsilon = abs_tol);
            let rel = if y.abs() > 1e-12 { (x - y).abs() / y.abs() } else { 0.0 };
            assert!(rel <= rel_tol, "rel error at {}: {} > {}", i, rel, rel_tol);
        }
    }

    /// Stage test: if C dump exists (run C with VLFEAT_SIFT_DUMP=1 first), compare.
    /// This is the incremental approach: octave base -> full octave -> dog -> grad.
    #[test]
    fn test_stage_octave_dump_if_present() {
        // Example: generate reference by running C code with env, then compare Rust output.
        // For now, just verify the loader works if file exists.
        if let Some(ref_data) = load_bin("/tmp/vl_sift_octave_first_base.bin") {
            // To compare, we would run the same input through Rust smooth/upsample and assert_approx_slices
            // Here we just check load succeeds and has reasonable size.
            assert!(!ref_data.is_empty());
            println!("Loaded C dump with {} floats", ref_data.len());
            // In full stage test: compute Rust equivalent, then assert_approx_slices(&rust_out, &ref_data, 1e-5, 1e-4);
        }
    }
}

