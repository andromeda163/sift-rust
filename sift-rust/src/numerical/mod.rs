//! Numerical building blocks ported from VLFeat SIFT (sift.c).
//! Each submodule corresponds to one verified primitive.

pub mod fast_expn;
pub mod smooth;
pub mod upsample;
pub mod downsample;

// Re-export the public API at the numerical:: level for convenience
pub use fast_expn::{fast_expn, fast_expn_init};
pub use smooth::{sift_smooth, vl_sift_smooth};
pub use upsample::copy_and_upsample_rows;
pub use downsample::copy_and_downsample;
