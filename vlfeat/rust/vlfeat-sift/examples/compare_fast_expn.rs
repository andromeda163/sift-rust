//! Example binary to generate Rust fast_expn values for comparison against C.
//! Run with: cargo run --example compare_fast_expn --release

use vlfeat_sift::fast_expn;

fn main() {
    // Test points: [0, 25] inclusive, plus some beyond
    let mut xs: Vec<f64> = (0..=250).map(|i| i as f64 * 0.1).collect(); // 0.0 .. 25.0 step 0.1
    xs.push(25.0); // exact boundary
    xs.push(25.0001);
    xs.push(30.0);

    println!("x,rust");
    for &x in &xs {
        let val = fast_expn(x);
        println!("{:.6},{:.12}", x, val);
    }
}