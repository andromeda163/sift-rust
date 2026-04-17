#!/bin/bash
set -euo pipefail

# test_against_c.sh
# Tests the Rust fast_expn (and can be extended to other core functions)
# against the equivalent C implementation for numerical accuracy.
#
# Usage: ./test_against_c.sh
# It builds the Rust example, compiles a minimal C driver (exact logic from vl/sift.c),
# runs both on the same grid of x values (including x==25.0 boundary), and compares
# with abs+rel tolerance.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Building Rust example ==="
cargo build --example compare_fast_expn --release 2>&1 | tail -5

echo "=== Writing and compiling C reference driver ==="
cat > /tmp/fast_expn_ref.c << 'EOF'
/* Minimal C driver for fast_expn reference (exact logic from vl/sift.c, post-fix) */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EXPN_SZ  256
#define EXPN_MAX 25.0
double expn_tab [EXPN_SZ + 1];

static void fast_expn_init(void) {
    for (int k = 0; k <= EXPN_SZ; ++k) {
        expn_tab[k] = exp( - (double)k * (EXPN_MAX / EXPN_SZ) );
    }
}

static double fast_expn(double x) {
    if (x >= EXPN_MAX) return 0.0;
    double xs = x * (EXPN_SZ / EXPN_MAX);
    int i = (int)floor(xs);
    double r = xs - i;
    double a = expn_tab[i];
    double b = (i < EXPN_SZ) ? expn_tab[i + 1] : a;
    return a + r * (b - a);
}

int main(void) {
    fast_expn_init();
    /* Same grid as Rust example */
    for (int i = 0; i <= 250; ++i) {
        double x = i * 0.1;
        printf("%.6f,%.12f\n", x, fast_expn(x));
    }
    printf("%.6f,%.12f\n", 25.0, fast_expn(25.0));
    printf("%.6f,%.12f\n", 25.0001, fast_expn(25.0001));
    printf("%.6f,%.12f\n", 30.0, fast_expn(30.0));
    return 0;
}
EOF

gcc -O2 -Wall -Wextra -o /tmp/fast_expn_ref /tmp/fast_expn_ref.c -lm
echo "C driver compiled."

echo "=== Running C reference ==="
/tmp/fast_expn_ref > /tmp/c_values.csv
echo "C values: $(wc -l < /tmp/c_values.csv) lines"

echo "=== Running Rust implementation ==="
cargo run --example compare_fast_expn --release 2>/dev/null | tail -n +2 > /tmp/rust_values.csv
echo "Rust values: $(wc -l < /tmp/rust_values.csv) lines"

echo "=== Comparing with tolerance (abs<=1e-9, rel<=1e-9) ==="
python3 - << 'PYEOF'
import csv, sys, math

def read_csv(path):
    xs, vs = [], []
    with open(path) as f:
        for line in f:
            if ',' not in line: continue
            x, v = line.strip().split(',')
            xs.append(float(x))
            vs.append(float(v))
    return xs, vs

cx, cv = read_csv('/tmp/c_values.csv')
rx, rv = read_csv('/tmp/rust_values.csv')

if len(cx) != len(rx):
    print(f"Length mismatch: C={len(cx)} Rust={len(rx)}")
    sys.exit(1)

abs_tol = 1e-9
rel_tol = 1e-9
max_abs_err = 0.0
max_rel_err = 0.0
failures = 0

for i, (x, c, r) in enumerate(zip(cx, cv, rv)):
    aerr = abs(c - r)
    rerr = aerr / abs(c) if abs(c) > 1e-12 else 0.0
    max_abs_err = max(max_abs_err, aerr)
    max_rel_err = max(max_rel_err, rerr)
    if aerr > abs_tol or rerr > rel_tol:
        failures += 1
        if failures <= 5:
            print(f"FAIL at x={x}: C={c:.12f} R={r:.12f} abs={aerr:.2e} rel={rerr:.2e}")

print(f"Max abs err: {max_abs_err:.2e}")
print(f"Max rel err: {max_rel_err:.2e}")
if failures == 0:
    print("PASS: Rust matches C within tolerance for all points (including x==25.0 boundary).")
    sys.exit(0)
else:
    print(f"FAIL: {failures} points exceeded tolerance.")
    sys.exit(1)
PYEOF

echo "=== Test script finished ==="