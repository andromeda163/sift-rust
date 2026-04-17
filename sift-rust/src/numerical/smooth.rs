//! On-the-fly Gaussian kernel + separable smoothing (replicate pad).
//! Matches C _vl_sift_smooth numerically within ~1e-6 (add order may differ slightly).

/// On-the-fly Gaussian kernel + separable smoothing (replicate pad).
/// Matches C _vl_sift_smooth numerically within ~1e-6 (add order may differ slightly).
pub fn vl_sift_smooth(
    output: &mut [f32],
    temp: &mut [f32],
    input: &[f32],
    width: usize,
    height: usize,
    sigma: f64,
) {
    let gfw = ((4.0 * sigma).ceil() as usize).max(1);
    let mut kernel = vec![0f32; 2 * gfw + 1];
    let mut acc = 0f32;
    for j in 0..=2 * gfw {
        let d = (j as isize - gfw as isize) as f32 / sigma as f32;
        kernel[j] = (-0.5 * d * d).exp() as f32;
        acc += kernel[j];
    }
    for v in kernel.iter_mut() {
        *v /= acc;
    }

    if gfw == 0 {
        output[..width * height].copy_from_slice(&input[..width * height]);
        return;
    }

    // Vertical pass (replicate pad) -> temp (row-major)
    for y in 0..height {
        for x in 0..width {
            let mut sum = 0f32;
            for (k, &w) in kernel.iter().enumerate() {
                let yy = (y as isize + (k as isize - gfw as isize)).clamp(0, height as isize - 1) as usize;
                sum += input[yy * width + x] * w;
            }
            temp[y * width + x] = sum;
        }
    }

    // Horizontal pass -> output
    for y in 0..height {
        for x in 0..width {
            let mut sum = 0f32;
            for (k, &w) in kernel.iter().enumerate() {
                let xx = (x as isize + (k as isize - gfw as isize)).clamp(0, width as isize - 1) as usize;
                sum += temp[y * width + xx] * w;
            }
            output[y * width + x] = sum;
        }
    }
}

/// Convenience wrapper.
pub fn sift_smooth(output: &mut [f32], temp: &mut [f32], input: &[f32], width: usize, height: usize, sigma: f64) {
    vl_sift_smooth(output, temp, input, width, height, sigma);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smooth_gaussian() {
        let w = 8usize;
        let h = 8usize;
        let mut input = vec![0f32; w*h];
        input[3 + 3*w] = 1.0;
        let mut temp = vec![0f32; w*h];
        let mut output = vec![0f32; w*h];
        vl_sift_smooth(&mut output, &mut temp, &input, w, h, 1.0);
        let sum: f32 = output.iter().sum();
        assert!((sum - 1.0).abs() < 0.01);
        assert!(*output.iter().max_by(|a,b| a.partial_cmp(b).unwrap()).unwrap() < 0.3);
    }
}
