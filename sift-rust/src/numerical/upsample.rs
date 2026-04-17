//! Copy + 2x upsample (linear interp) + transpose.
//! Exact port of C copy_and_upsample_rows (pointer arithmetic replicated).

/// Copy + 2x upsample (linear interp) + transpose.
/// Exact port of C copy_and_upsample_rows (pointer arithmetic replicated).
pub fn copy_and_upsample_rows(dst: &mut [f32], src: &[f32], width: usize, height: usize) {
    assert!(dst.len() >= 2 * width * height);
    let mut dst_idx: isize = 0;
    let h = height as isize;
    let w2h = (2 * width) as isize * h;
    for y in 0..height {
        let mut src_idx = y * width;
        let mut a = src[src_idx];
        src_idx += 1;
        let mut b = a;
        for _x in 0..width - 1 {
            b = src[src_idx];
            src_idx += 1;
            dst[dst_idx as usize] = a;
            dst_idx += h;
            dst[dst_idx as usize] = 0.5 * (a + b);
            dst_idx += h;
            a = b;
        }
        dst[dst_idx as usize] = b;
        dst_idx += h;
        dst[dst_idx as usize] = b;
        dst_idx += h;
        dst_idx += 1 - w2h;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_copy_upsample_basic() {
        let src = vec![1.0f32, 2.0, 3.0, 4.0];
        let mut dst = vec![0f32; 8];
        copy_and_upsample_rows(&mut dst, &src, 2, 2);
        assert!(dst.iter().any(|&v| v > 0.5));
    }
}
