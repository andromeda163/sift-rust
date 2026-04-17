//! Downsample by 2^d (nearest, no average).
//! Exact port of C copy_and_downsample.

/// Downsample by 2^d (nearest, no average).
/// Exact port of C copy_and_downsample.
pub fn copy_and_downsample(dst: &mut [f32], src: &[f32], width: usize, height: usize, d: usize) {
    let step = 1usize << d;
    let mut di = 0usize;
    for y in (0..height).step_by(step) {
        let srcrowp = y * width;
        for x in (0..width).step_by(step) {
            if x < width - (step - 1) {
                dst[di] = src[srcrowp + x];
                di += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_downsample() {
        let src = vec![1f32,2.,3.,4., 5.,6.,7.,8.];
        let mut dst = vec![0f32; 2];
        copy_and_downsample(&mut dst, &src, 4, 2, 1);
        assert_eq!(dst[0], 1.0);
    }
}
