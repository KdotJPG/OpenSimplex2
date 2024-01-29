/*!
    K.jpg's OpenSimplex 2, smooth variant ("SuperSimplex")
*/

use std::{num::Wrapping, sync::Once};

const PRIME_X: i64 = 0x5205402B9270C86F;
const PRIME_Y: i64 = 0x598CD327003817B5;
const PRIME_Z: i64 = 0x5BCC226E9FA0BACB;
const PRIME_W: i64 = 0x56CC5227E58F554B;
const HASH_MULTIPLIER: i64 = 0x53A3F72DEEC546F5;
const SEED_FLIP_3D: i64 = -0x52D547B2E96ED629;

const ROOT2OVER2: f64 = 0.7071067811865476;
const SKEW_2D: f64 = 0.366025403784439;
const UNSKEW_2D: f64 = -0.21132486540518713;

const ROOT3OVER3: f64 = 0.577350269189626;
const FALLBACK_ROTATE3: f64 = 2.0 / 3.0;
const ROTATE3_ORTHOGONALIZER: f64 = UNSKEW_2D;

const SKEW_4D: f32 = 0.309016994374947;
const UNSKEW_4D: f32 = -0.138196601125011;

const N_GRADS_2D_EXPONENT: i32 = 7;
const N_GRADS_3D_EXPONENT: i32 = 8;
const N_GRADS_4D_EXPONENT: i32 = 9;
const N_GRADS_2D: i32 = 1 << N_GRADS_2D_EXPONENT;
const N_GRADS_3D: i32 = 1 << N_GRADS_3D_EXPONENT;
const N_GRADS_4D: i32 = 1 << N_GRADS_4D_EXPONENT;

const NORMALIZER_2D: f64 = 0.05481866495625118;
const NORMALIZER_3D: f64 = 0.2781926117527186;
const NORMALIZER_4D: f64 = 0.11127401889945551;

const RSQUARED_2D: f32 = 2.0 / 3.0;
const RSQUARED_3D: f32 = 3.0 / 4.0;
const RSQUARED_4D: f32 = 4.0 / 5.0;

/*
    Noise Evaluators
*/

/**
    2D OpenSimplex2S/SuperSimplex noise, standard lattice orientation.
*/
pub fn noise2(seed: i64, x: f64, y: f64) -> f32 {
    // Get points for A2* lattice
    let s = SKEW_2D * (x + y);
    let xs = x + s;
    let ys = y + s;

    noise2_UnskewedBase(seed, xs, ys)
}

/**
    2D OpenSimplex2S/SuperSimplex noise, with Y pointing down the main diagonal.
    Might be better for a 2D sandbox style game, where Y is vertical.
    Probably slightly less optimal for heightmaps or continent maps,
    unless your map is centered around an equator. It's a slight
    difference, but the option is here to make it easy.
*/
pub fn noise2_ImproveX(seed: i64, x: f64, y: f64) -> f32 {
    // Skew transform and rotation baked into one.
    let xx = x * ROOT2OVER2;
    let yy = y * (ROOT2OVER2 * (1.0 + 2.0 * SKEW_2D));

    noise2_UnskewedBase(seed, yy + xx, yy - xx)
}

/**
    2D  OpenSimplex2S/SuperSimplex noise base.
*/
fn noise2_UnskewedBase(seed: i64, xs: f64, ys: f64) -> f32 {
    let seed = Wrapping(seed);

    // Get base points and offsets.
    let xsb = fastFloor(xs);
    let ysb = fastFloor(ys);
    let xi = (xs - xsb as f64) as f32;
    let yi = (ys - ysb as f64) as f32;

    // Prime pre-multiplication for hash.
    let xsbp = Wrapping(xsb as i64) * Wrapping(PRIME_X);
    let ysbp = Wrapping(ysb as i64) * Wrapping(PRIME_Y);

    // Unskew.
    let t = (xi + yi) * UNSKEW_2D as f32;
    let dx0 = xi + t;
    let dy0 = yi + t;

    // First vertex.
    let a0 = RSQUARED_2D - dx0 * dx0 - dy0 * dy0;
    let mut value = (a0 * a0) * (a0 * a0) * grad2(seed, xsbp, ysbp, dx0, dy0);

    // Second vertex.
    let a1 = (2.0 * (1.0 + 2.0 * UNSKEW_2D) * (1.0 / UNSKEW_2D + 2.0)) as f32 * t
        + ((-2.0 * (1.0 + 2.0 * UNSKEW_2D) * (1.0 + 2.0 * UNSKEW_2D)) as f32 + a0);
    let dx1 = dx0 - (1.0 + 2.0 * UNSKEW_2D) as f32;
    let dy1 = dy0 - (1.0 + 2.0 * UNSKEW_2D) as f32;
    value += (a1 * a1)
        * (a1 * a1)
        * grad2(
            seed,
            xsbp + Wrapping(PRIME_X),
            ysbp + Wrapping(PRIME_Y),
            dx1,
            dy1,
        );

    // Third and fourth vertices.
    // Nested conditionals were faster than compact bit logic/arithmetic.
    let xmyi = xi - yi;
    if t < UNSKEW_2D as f32 {
        if xi + xmyi > 1.0 {
            let dx2 = dx0 - (3.0 * UNSKEW_2D + 2.0) as f32;
            let dy2 = dy0 - (3.0 * UNSKEW_2D + 1.0) as f32;
            let a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if a2 > 0.0 {
                value += (a2 * a2)
                    * (a2 * a2)
                    * grad2(
                        seed,
                        xsbp + Wrapping(PRIME_X << 1),
                        ysbp + Wrapping(PRIME_Y),
                        dx2,
                        dy2,
                    );
            }
        } else {
            let dx2 = dx0 - UNSKEW_2D as f32;
            let dy2 = dy0 - (UNSKEW_2D + 1.0) as f32;
            let a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if a2 > 0.0 {
                value +=
                    (a2 * a2) * (a2 * a2) * grad2(seed, xsbp, ysbp + Wrapping(PRIME_Y), dx2, dy2);
            }
        }

        if yi - xmyi > 1.0 {
            let dx3 = dx0 - (3.0 * UNSKEW_2D + 1.0) as f32;
            let dy3 = dy0 - (3.0 * UNSKEW_2D + 2.0) as f32;
            let a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
            if a3 > 0.0 {
                value += (a3 * a3)
                    * (a3 * a3)
                    * grad2(
                        seed,
                        xsbp + Wrapping(PRIME_X),
                        ysbp + Wrapping(PRIME_Y << 1),
                        dx3,
                        dy3,
                    );
            }
        } else {
            let dx3 = dx0 - (UNSKEW_2D + 1.0) as f32;
            let dy3 = dy0 - UNSKEW_2D as f32;
            let a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
            if a3 > 0.0 {
                value +=
                    (a3 * a3) * (a3 * a3) * grad2(seed, xsbp + Wrapping(PRIME_X), ysbp, dx3, dy3);
            }
        }
    } else {
        if xi + xmyi < 0.0 {
            let dx2 = dx0 + (1.0 + UNSKEW_2D) as f32;
            let dy2 = dy0 + UNSKEW_2D as f32;
            let a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if a2 > 0.0 {
                value +=
                    (a2 * a2) * (a2 * a2) * grad2(seed, xsbp - Wrapping(PRIME_X), ysbp, dx2, dy2);
            }
        } else {
            let dx2 = dx0 - (UNSKEW_2D + 1.0) as f32;
            let dy2 = dy0 - UNSKEW_2D as f32;
            let a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if a2 > 0.0 {
                value +=
                    (a2 * a2) * (a2 * a2) * grad2(seed, xsbp + Wrapping(PRIME_X), ysbp, dx2, dy2);
            }
        }

        if yi < xmyi {
            let dx2 = dx0 + UNSKEW_2D as f32;
            let dy2 = dy0 + (UNSKEW_2D + 1.0) as f32;
            let a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if a2 > 0.0 {
                value +=
                    (a2 * a2) * (a2 * a2) * grad2(seed, xsbp, ysbp - Wrapping(PRIME_Y), dx2, dy2);
            }
        } else {
            let dx2 = dx0 - UNSKEW_2D as f32;
            let dy2 = dy0 - (UNSKEW_2D + 1.0) as f32;
            let a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if a2 > 0.0 {
                value +=
                    (a2 * a2) * (a2 * a2) * grad2(seed, xsbp, ysbp + Wrapping(PRIME_Y), dx2, dy2);
            }
        }
    }

    value
}

/**
    3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Y).
    Recommended for 3D terrain and time-varied animations.
    The Z coordinate should always be the "different" coordinate in whatever your use case is.
    If Y is vertical in world coordinates, call noise3_ImproveXZ(x, z, Y) or use noise3_XZBeforeY.
    If Z is vertical in world coordinates, call noise3_ImproveXZ(x, y, Z).
    For a time varied animation, call noise3_ImproveXY(x, y, T).
*/
pub fn noise3_ImproveXY(seed: i64, x: f64, y: f64, z: f64) -> f32 {
    // Re-orient the cubic lattices without skewing, so Z points up the main lattice diagonal,
    // and the planes formed by XY are moved far out of alignment with the cube faces.
    // Orthonormal rotation. Not a skew transform.
    let xy = x + y;
    let s2 = xy * ROTATE3_ORTHOGONALIZER;
    let zz = z * ROOT3OVER3;
    let xr = x + s2 + zz;
    let yr = y + s2 + zz;
    let zr = xy * -ROOT3OVER3 + zz;

    // Evaluate both lattices to form a BCC lattice.
    noise3_UnrotatedBase(seed, xr, yr, zr)
}

/**
    3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Z).
    Recommended for 3D terrain and time-varied animations.
    The Y coordinate should always be the "different" coordinate in whatever your use case is.
    If Y is vertical in world coordinates, call noise3_ImproveXZ(x, Y, z).
    If Z is vertical in world coordinates, call noise3_ImproveXZ(x, Z, y) or use noise3_ImproveXY.
    For a time varied animation, call noise3_ImproveXZ(x, T, y) or use noise3_ImproveXY.
*/
pub fn noise3_ImproveXZ(seed: i64, x: f64, y: f64, z: f64) -> f32 {
    // Re-orient the cubic lattices without skewing, so Y points up the main lattice diagonal,
    // and the planes formed by XZ are moved far out of alignment with the cube faces.
    // Orthonormal rotation. Not a skew transform.
    let xz = x + z;
    let s2 = xz * -0.211324865405187;
    let yy = y * ROOT3OVER3;
    let xr = x + s2 + yy;
    let zr = z + s2 + yy;
    let yr = xz * -ROOT3OVER3 + yy;

    // Evaluate both lattices to form a BCC lattice.
    noise3_UnrotatedBase(seed, xr, yr, zr)
}

/**
    3D OpenSimplex2S/SuperSimplex noise, fallback rotation option
    Use noise3_ImproveXY or noise3_ImproveXZ instead, wherever appropriate.
    They have less diagonal bias. This function's best use is as a fallback.
*/
pub fn noise3_Fallback(seed: i64, x: f64, y: f64, z: f64) -> f32 {
    // Re-orient the cubic lattices via rotation, to produce a familiar look.
    // Orthonormal rotation. Not a skew transform.
    let r = FALLBACK_ROTATE3 * (x + y + z);
    let xr = r - x;
    let yr = r - y;
    let zr = r - z;

    // Evaluate both lattices to form a BCC lattice.
    noise3_UnrotatedBase(seed, xr, yr, zr)
}

/**
    Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
    Lookup table implementation inspired by DigitalShadow.
    It was actually faster to narrow down the points in the loop itself,
    than to build up the index with enough info to isolate 8 points.
*/
fn noise3_UnrotatedBase(seed: i64, xr: f64, yr: f64, zr: f64) -> f32 {
    let seed = Wrapping(seed);

    // Get base points and offsets.
    let xrb = fastFloor(xr);
    let yrb = fastFloor(yr);
    let zrb = fastFloor(zr);
    let xi = (xr - xrb as f64) as f32;
    let yi = (yr - yrb as f64) as f32;
    let zi = (zr - zrb as f64) as f32;

    // Prime pre-multiplication for hash. Also flip seed for second lattice copy.
    let xrbp = Wrapping(xrb as i64) * Wrapping(PRIME_X);
    let yrbp = Wrapping(yrb as i64) * Wrapping(PRIME_Y);
    let zrbp = Wrapping(zrb as i64) * Wrapping(PRIME_Z);
    let seed2 = seed ^ Wrapping(SEED_FLIP_3D);

    // -1 if positive, 0 if negative.
    let xNMask = (-0.5 - xi) as i32;
    let yNMask = (-0.5 - yi) as i32;
    let zNMask = (-0.5 - zi) as i32;

    // First vertex.
    let x0 = xi + xNMask as f32;
    let y0 = yi + yNMask as f32;
    let z0 = zi + zNMask as f32;
    let a0 = RSQUARED_3D - x0 * x0 - y0 * y0 - z0 * z0;
    let mut value = (a0 * a0)
        * (a0 * a0)
        * grad3(
            seed,
            xrbp + (Wrapping(xNMask as i64) & Wrapping(PRIME_X)),
            yrbp + (Wrapping(yNMask as i64) & Wrapping(PRIME_Y)),
            zrbp + (Wrapping(zNMask as i64) & Wrapping(PRIME_Z)),
            x0,
            y0,
            z0,
        );

    // Second vertex.
    let x1 = xi - 0.5;
    let y1 = yi - 0.5;
    let z1 = zi - 0.5;
    let a1 = RSQUARED_3D - x1 * x1 - y1 * y1 - z1 * z1;
    value += (a1 * a1)
        * (a1 * a1)
        * grad3(
            seed2,
            xrbp + Wrapping(PRIME_X),
            yrbp + Wrapping(PRIME_Y),
            zrbp + Wrapping(PRIME_Z),
            x1,
            y1,
            z1,
        );

    // Shortcuts for building the remaining falloffs.
    // Derived by subtracting the polynomials with the offsets plugged in.
    let xAFlipMask0 = ((xNMask | 1) << 1) as f32 * x1;
    let yAFlipMask0 = ((yNMask | 1) << 1) as f32 * y1;
    let zAFlipMask0 = ((zNMask | 1) << 1) as f32 * z1;
    let xAFlipMask1 = (-2 - (xNMask << 2)) as f32 * x1 - 1.0;
    let yAFlipMask1 = (-2 - (yNMask << 2)) as f32 * y1 - 1.0;
    let zAFlipMask1 = (-2 - (zNMask << 2)) as f32 * z1 - 1.0;

    let mut skip5 = false;
    let a2 = xAFlipMask0 + a0;
    if a2 > 0.0 {
        let x2 = x0 - (xNMask | 1) as f32;
        let y2 = y0;
        let z2 = z0;
        value += (a2 * a2)
            * (a2 * a2)
            * grad3(
                seed,
                xrbp + (Wrapping(!xNMask as i64) & Wrapping(PRIME_X)),
                yrbp + (Wrapping(yNMask as i64) & Wrapping(PRIME_Y)),
                zrbp + (Wrapping(zNMask as i64) & Wrapping(PRIME_Z)),
                x2,
                y2,
                z2,
            );
    } else {
        let a3 = yAFlipMask0 + zAFlipMask0 + a0;
        if a3 > 0.0 {
            let x3 = x0;
            let y3 = y0 - (yNMask | 1) as f32;
            let z3 = z0 - (zNMask | 1) as f32;
            value += (a3 * a3)
                * (a3 * a3)
                * grad3(
                    seed,
                    xrbp + (Wrapping(xNMask as i64) & Wrapping(PRIME_X)),
                    yrbp + (Wrapping(!yNMask as i64) & Wrapping(PRIME_Y)),
                    zrbp + (Wrapping(!zNMask as i64) & Wrapping(PRIME_Z)),
                    x3,
                    y3,
                    z3,
                );
        }

        let a4 = xAFlipMask1 + a1;
        if a4 > 0.0 {
            let x4 = (xNMask | 1) as f32 + x1;
            let y4 = y1;
            let z4 = z1;
            value += (a4 * a4)
                * (a4 * a4)
                * grad3(
                    seed2,
                    xrbp + (Wrapping(xNMask as i64) & (Wrapping(PRIME_X) << 1)),
                    yrbp + Wrapping(PRIME_Y),
                    zrbp + Wrapping(PRIME_Z),
                    x4,
                    y4,
                    z4,
                );
            skip5 = true;
        }
    }

    let mut skip9 = false;
    let a6 = yAFlipMask0 + a0;
    if a6 > 0.0 {
        let x6 = x0;
        let y6 = y0 - (yNMask | 1) as f32;
        let z6 = z0;
        value += (a6 * a6)
            * (a6 * a6)
            * grad3(
                seed,
                xrbp + (Wrapping(xNMask as i64) & Wrapping(PRIME_X)),
                yrbp + (Wrapping(!yNMask as i64) & Wrapping(PRIME_Y)),
                zrbp + (Wrapping(zNMask as i64) & Wrapping(PRIME_Z)),
                x6,
                y6,
                z6,
            );
    } else {
        let a7 = xAFlipMask0 + zAFlipMask0 + a0;
        if a7 > 0.0 {
            let x7 = x0 - (xNMask | 1) as f32;
            let y7 = y0;
            let z7 = z0 - (zNMask | 1) as f32;
            value += (a7 * a7)
                * (a7 * a7)
                * grad3(
                    seed,
                    xrbp + (Wrapping(!xNMask as i64) & Wrapping(PRIME_X)),
                    yrbp + (Wrapping(yNMask as i64) & Wrapping(PRIME_Y)),
                    zrbp + (Wrapping(!zNMask as i64) & Wrapping(PRIME_Z)),
                    x7,
                    y7,
                    z7,
                );
        }

        let a8 = yAFlipMask1 + a1;
        if a8 > 0.0 {
            let x8 = x1;
            let y8 = (yNMask | 1) as f32 + y1;
            let z8 = z1;
            value += (a8 * a8)
                * (a8 * a8)
                * grad3(
                    seed2,
                    xrbp + Wrapping(PRIME_X),
                    yrbp + (Wrapping(yNMask as i64) & (Wrapping(PRIME_Y) << 1)),
                    zrbp + Wrapping(PRIME_Z),
                    x8,
                    y8,
                    z8,
                );
            skip9 = true;
        }
    }

    let mut skipD = false;
    let aA = zAFlipMask0 + a0;
    if aA > 0.0 {
        let xA = x0;
        let yA = y0;
        let zA = z0 - (zNMask | 1) as f32;
        value += (aA * aA)
            * (aA * aA)
            * grad3(
                seed,
                xrbp + (Wrapping(xNMask as i64) & Wrapping(PRIME_X)),
                yrbp + (Wrapping(yNMask as i64) & Wrapping(PRIME_Y)),
                zrbp + (Wrapping(!zNMask as i64) & Wrapping(PRIME_Z)),
                xA,
                yA,
                zA,
            );
    } else {
        let aB = xAFlipMask0 + yAFlipMask0 + a0;
        if aB > 0.0 {
            let xB = x0 - (xNMask | 1) as f32;
            let yB = y0 - (yNMask | 1) as f32;
            let zB = z0;
            value += (aB * aB)
                * (aB * aB)
                * grad3(
                    seed,
                    xrbp + (Wrapping(!xNMask as i64) & Wrapping(PRIME_X)),
                    yrbp + (Wrapping(!yNMask as i64) & Wrapping(PRIME_Y)),
                    zrbp + (Wrapping(zNMask as i64) & Wrapping(PRIME_Z)),
                    xB,
                    yB,
                    zB,
                );
        }

        let aC = zAFlipMask1 + a1;
        if aC > 0.0 {
            let xC = x1;
            let yC = y1;
            let zC = (zNMask | 1) as f32 + z1;
            value += (aC * aC)
                * (aC * aC)
                * grad3(
                    seed2,
                    xrbp + Wrapping(PRIME_X),
                    yrbp + Wrapping(PRIME_Y),
                    zrbp + (Wrapping(zNMask as i64) & (Wrapping(PRIME_Z) << 1)),
                    xC,
                    yC,
                    zC,
                );
            skipD = true;
        }
    }

    if !skip5 {
        let a5 = yAFlipMask1 + zAFlipMask1 + a1;
        if a5 > 0.0 {
            let x5 = x1;
            let y5 = (yNMask | 1) as f32 + y1;
            let z5 = (zNMask | 1) as f32 + z1;
            value += (a5 * a5)
                * (a5 * a5)
                * grad3(
                    seed2,
                    xrbp + Wrapping(PRIME_X),
                    yrbp + (Wrapping(yNMask as i64) & (Wrapping(PRIME_Y) << 1)),
                    zrbp + (Wrapping(zNMask as i64) & (Wrapping(PRIME_Z) << 1)),
                    x5,
                    y5,
                    z5,
                );
        }
    }

    if !skip9 {
        let a9 = xAFlipMask1 + zAFlipMask1 + a1;
        if a9 > 0.0 {
            let x9 = (xNMask | 1) as f32 + x1;
            let y9 = y1;
            let z9 = (zNMask | 1) as f32 + z1;
            value += (a9 * a9)
                * (a9 * a9)
                * grad3(
                    seed2,
                    xrbp + (Wrapping(xNMask as i64) & (Wrapping(PRIME_X) << 1)),
                    yrbp + Wrapping(PRIME_Y),
                    zrbp + (Wrapping(zNMask as i64) & (Wrapping(PRIME_Z) << 1)),
                    x9,
                    y9,
                    z9,
                );
        }
    }

    if !skipD {
        let aD = xAFlipMask1 + yAFlipMask1 + a1;
        if aD > 0.0 {
            let xD = (xNMask | 1) as f32 + x1;
            let yD = (yNMask | 1) as f32 + y1;
            let zD = z1;
            value += (aD * aD)
                * (aD * aD)
                * grad3(
                    seed2,
                    xrbp + (Wrapping(xNMask as i64) & (Wrapping(PRIME_X) << 1)),
                    yrbp + (Wrapping(yNMask as i64) & (Wrapping(PRIME_Y) << 1)),
                    zrbp + Wrapping(PRIME_Z),
                    xD,
                    yD,
                    zD,
                );
        }
    }

    value
}

/**
    4D SuperSimplex noise, with XYZ oriented like noise3_ImproveXY
    and W for an extra degree of freedom. W repeats eventually.
    Recommended for time-varied animations which texture a 3D object (W=time)
    in a space where Z is vertical
*/
pub fn noise4_ImproveXYZ_ImproveXY(seed: i64, x: f64, y: f64, z: f64, w: f64) -> f32 {
    let xy = x + y;
    let s2 = xy * -0.21132486540518699998;
    let zz = z * 0.28867513459481294226;
    let ww = w * 1.118033988749894;
    let xr = x + (zz + ww + s2);
    let yr = y + (zz + ww + s2);
    let zr = xy * -0.57735026918962599998 + (zz + ww);
    let wr = z * -0.866025403784439 + ww;

    noise4_UnskewedBase(seed, xr, yr, zr, wr)
}

/**
    4D SuperSimplex noise, with XYZ oriented like noise3_ImproveXZ
    and W for an extra degree of freedom. W repeats eventually.
    Recommended for time-varied animations which texture a 3D object (W=time)
    in a space where Y is vertical
*/
pub fn noise4_ImproveXYZ_ImproveXZ(seed: i64, x: f64, y: f64, z: f64, w: f64) -> f32 {
    let xz = x + z;
    let s2 = xz * -0.21132486540518699998;
    let yy = y * 0.28867513459481294226;
    let ww = w * 1.118033988749894;
    let xr = x + (yy + ww + s2);
    let zr = z + (yy + ww + s2);
    let yr = xz * -0.57735026918962599998 + (yy + ww);
    let wr = y * -0.866025403784439 + ww;

    noise4_UnskewedBase(seed, xr, yr, zr, wr)
}

/**
    4D SuperSimplex noise, with XYZ oriented like noise3_Fallback
    and W for an extra degree of freedom. W repeats eventually.
    Recommended for time-varied animations which texture a 3D object (W=time)
    where there isn't a clear distinction between horizontal and vertical
*/
pub fn noise4_ImproveXYZ(seed: i64, x: f64, y: f64, z: f64, w: f64) -> f32 {
    let xyz = x + y + z;
    let ww = w * 1.118033988749894;
    let s2 = xyz * -0.16666666666666666 + ww;
    let xs = x + s2;
    let ys = y + s2;
    let zs = z + s2;
    let ws = -0.5 * xyz + ww;

    noise4_UnskewedBase(seed, xs, ys, zs, ws)
}

/**
    4D SuperSimplex noise, with XY and ZW forming orthogonal triangular-based planes.
    Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
    Recommended for noise(x, y, sin(time), cos(time)) trick.
*/
pub fn noise4_ImproveXY_ImproveZW(seed: i64, x: f64, y: f64, z: f64, w: f64) -> f32 {
    let s2 = (x + y) * -0.28522513987434876941 + (z + w) * 0.83897065470611435718;
    let t2 = (z + w) * 0.21939749883706435719 + (x + y) * -0.48214856493302476942;
    let xs = x + s2;
    let ys = y + s2;
    let zs = z + t2;
    let ws = w + t2;

    noise4_UnskewedBase(seed, xs, ys, zs, ws)
}

/**
    4D SuperSimplex noise, fallback lattice orientation.
*/
pub fn noise4_Fallback(seed: i64, x: f64, y: f64, z: f64, w: f64) -> f32 {
    // Get points for A4 lattice
    let s = SKEW_4D as f64 * (x + y + z + w);
    let xs = x + s;
    let ys = y + s;
    let zs = z + s;
    let ws = w + s;

    noise4_UnskewedBase(seed, xs, ys, zs, ws)
}

/**
    4D SuperSimplex noise base.
    Using ultra-simple 4x4x4x4 lookup partitioning.
    This isn't as elegant or SIMD/GPU/etc. portable as other approaches,
    but it competes performance-wise with optimized 2014 OpenSimplex.
*/
fn noise4_UnskewedBase(seed: i64, xs: f64, ys: f64, zs: f64, ws: f64) -> f32 {
    let seed = Wrapping(seed);

    // Get base points and offsets
    let xsb = fastFloor(xs);
    let ysb = fastFloor(ys);
    let zsb = fastFloor(zs);
    let wsb = fastFloor(ws);
    let xsi = (xs - xsb as f64) as f32;
    let ysi = (ys - ysb as f64) as f32;
    let zsi = (zs - zsb as f64) as f32;
    let wsi = (ws - wsb as f64) as f32;

    // Unskewed offsets
    let ssi = (xsi + ysi + zsi + wsi) * UNSKEW_4D;
    let xi = xsi + ssi;
    let yi = ysi + ssi;
    let zi = zsi + ssi;
    let wi = wsi + ssi;

    // Prime pre-multiplication for hash.
    let xsvp = Wrapping(xsb as i64) * Wrapping(PRIME_X);
    let ysvp = Wrapping(ysb as i64) * Wrapping(PRIME_Y);
    let zsvp = Wrapping(zsb as i64) * Wrapping(PRIME_Z);
    let wsvp = Wrapping(wsb as i64) * Wrapping(PRIME_W);

    // Index into initial table.
    let index = ((fastFloor(xs * 4.0) & 3) << 0)
        | ((fastFloor(ys * 4.0) & 3) << 2)
        | ((fastFloor(zs * 4.0) & 3) << 4)
        | ((fastFloor(ws * 4.0) & 3) << 6);

    // Point contributions
    let staticData = getStaticData();
    let mut value = 0.0;
    let secondaryIndexStartAndStop = staticData.lookup4DA[index as usize];
    let secondaryIndexStart = secondaryIndexStartAndStop & 0xFFFF;
    let secondaryIndexStop = secondaryIndexStartAndStop >> 16;
    for i in secondaryIndexStart..secondaryIndexStop {
        let c = &staticData.lookup4DB[i];
        let dx = xi + c.dx;
        let dy = yi + c.dy;
        let dz = zi + c.dz;
        let dw = wi + c.dw;
        let mut a = (dx * dx + dy * dy) + (dz * dz + dw * dw);
        if a < RSQUARED_4D {
            a -= RSQUARED_4D;
            a *= a;
            value += a
                * a
                * grad4(
                    seed,
                    xsvp + Wrapping(c.xsvp),
                    ysvp + Wrapping(c.ysvp),
                    zsvp + Wrapping(c.zsvp),
                    wsvp + Wrapping(c.wsvp),
                    dx,
                    dy,
                    dz,
                    dw,
                );
        }
    }
    value
}

/*
    Utility
*/

fn grad2(seed: Wrapping<i64>, xsvp: Wrapping<i64>, ysvp: Wrapping<i64>, dx: f32, dy: f32) -> f32 {
    let mut hash = seed ^ xsvp ^ ysvp;
    hash *= HASH_MULTIPLIER;
    hash ^= hash.0 >> (64 - N_GRADS_2D_EXPONENT + 1);
    let gi = (hash.0 as i32 & ((N_GRADS_2D - 1) << 1)) as usize;
    let grads = &getStaticData().gradients2D;
    grads[gi | 0] * dx + grads[gi | 1] * dy
}

fn grad3(
    seed: Wrapping<i64>,
    xrvp: Wrapping<i64>,
    yrvp: Wrapping<i64>,
    zrvp: Wrapping<i64>,
    dx: f32,
    dy: f32,
    dz: f32,
) -> f32 {
    let mut hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
    hash *= HASH_MULTIPLIER;
    hash ^= hash.0 >> (64 - N_GRADS_3D_EXPONENT + 2);
    let gi = (hash.0 as i32 & ((N_GRADS_3D - 1) << 2)) as usize;
    let grads = &getStaticData().gradients3D;
    grads[gi | 0] * dx + grads[gi | 1] * dy + grads[gi | 2] * dz
}

fn grad4(
    seed: Wrapping<i64>,
    xsvp: Wrapping<i64>,
    ysvp: Wrapping<i64>,
    zsvp: Wrapping<i64>,
    wsvp: Wrapping<i64>,
    dx: f32,
    dy: f32,
    dz: f32,
    dw: f32,
) -> f32 {
    let mut hash = seed ^ (xsvp ^ ysvp) ^ (zsvp ^ wsvp);
    hash *= HASH_MULTIPLIER;
    hash ^= hash.0 >> (64 - N_GRADS_4D_EXPONENT + 2);
    let gi = (hash.0 as i32 & ((N_GRADS_4D - 1) << 2)) as usize;
    let grads = &getStaticData().gradients4D;
    (grads[gi | 0] * dx + grads[gi | 1] * dy) + (grads[gi | 2] * dz + grads[gi | 3] * dw)
}

fn fastFloor(x: f64) -> i32 {
    let xi = x as i32;
    if x < xi as f64 {
        xi - 1
    } else {
        xi
    }
}

/*
    Lookup Tables & Gradients
*/

#[derive(Clone, Default)]
struct LatticeVertex4D {
    pub dx: f32,
    pub dy: f32,
    pub dz: f32,
    pub dw: f32,
    pub xsvp: i64,
    pub ysvp: i64,
    pub zsvp: i64,
    pub wsvp: i64,
}

impl LatticeVertex4D {
    pub fn new(xsv: i32, ysv: i32, zsv: i32, wsv: i32) -> Self {
        let ssv = (xsv + ysv + zsv + wsv) as f32 * UNSKEW_4D;
        Self {
            xsvp: (Wrapping(xsv as i64) * Wrapping(PRIME_X)).0,
            ysvp: (Wrapping(ysv as i64) * Wrapping(PRIME_Y)).0,
            zsvp: (Wrapping(zsv as i64) * Wrapping(PRIME_Z)).0,
            wsvp: (Wrapping(wsv as i64) * Wrapping(PRIME_W)).0,
            dx: -xsv as f32 - ssv,
            dy: -ysv as f32 - ssv,
            dz: -zsv as f32 - ssv,
            dw: -wsv as f32 - ssv,
        }
    }
}

struct StaticData {
    gradients2D: Vec<f32>,
    gradients3D: Vec<f32>,
    gradients4D: Vec<f32>,

    lookup4DA: Vec<usize>,
    lookup4DB: Vec<LatticeVertex4D>,
}

static mut STATIC_DATA: (Once, Option<StaticData>) = (Once::new(), None);

fn getStaticData() -> &'static StaticData {
    unsafe {
        STATIC_DATA.0.call_once(|| {
            STATIC_DATA.1 = Some(initStaticData());
        });
        STATIC_DATA.1.as_ref().unwrap()
    }
}

fn initStaticData() -> StaticData {
    let gradients2D: Vec<_> = GRAD2_SRC
        .into_iter()
        .map(|v| (v / NORMALIZER_2D) as f32)
        .collect::<Vec<_>>() // cache divisions
        .into_iter()
        .cycle()
        .take((N_GRADS_2D * 2) as usize)
        .collect();

    let gradients3D: Vec<_> = GRAD3_SRC
        .into_iter()
        .map(|v| (v / NORMALIZER_3D) as f32)
        .collect::<Vec<_>>() // cache divisions
        .into_iter()
        .cycle()
        .take((N_GRADS_3D * 4) as usize)
        .collect();

    let gradients4D: Vec<_> = GRAD4_SRC
        .into_iter()
        .map(|v| (v / NORMALIZER_4D) as f32)
        .collect::<Vec<_>>() // cache divisions
        .into_iter()
        .cycle()
        .take((N_GRADS_4D * 4) as usize)
        .collect();

    let nLatticeVerticesTotal = LOOKUP_4D_VERTEX_CODES.iter().map(|v| v.len()).sum();
    let latticeVerticesByCode: Vec<_> = (0..256)
        .map(|i| {
            let cx = ((i >> 0) & 3) - 1;
            let cy = ((i >> 2) & 3) - 1;
            let cz = ((i >> 4) & 3) - 1;
            let cw = ((i >> 6) & 3) - 1;
            LatticeVertex4D::new(cx, cy, cz, cw)
        })
        .collect();
    let mut lookup4DA = vec![0; 256];
    let mut lookup4DB = vec![Default::default(); nLatticeVerticesTotal];
    let mut j = 0;
    for i in 0..256 {
        lookup4DA[i] = j | ((j + LOOKUP_4D_VERTEX_CODES[i].len()) << 16);
        for k in 0..LOOKUP_4D_VERTEX_CODES[i].len() {
            lookup4DB[j] = latticeVerticesByCode[LOOKUP_4D_VERTEX_CODES[i][k] as usize].clone();
            j += 1;
        }
    }

    StaticData {
        gradients2D,
        gradients3D,
        gradients4D,
        lookup4DA,
        lookup4DB,
    }
}

#[rustfmt::skip]
const GRAD2_SRC: &[f64] = &[
    0.38268343236509,   0.923879532511287,
    0.923879532511287,  0.38268343236509,
    0.923879532511287, -0.38268343236509,
    0.38268343236509,  -0.923879532511287,
    -0.38268343236509,  -0.923879532511287,
    -0.923879532511287, -0.38268343236509,
    -0.923879532511287,  0.38268343236509,
    -0.38268343236509,   0.923879532511287,
    //-------------------------------------//
    0.130526192220052,  0.99144486137381,
    0.608761429008721,  0.793353340291235,
    0.793353340291235,  0.608761429008721,
    0.99144486137381,   0.130526192220051,
    0.99144486137381,  -0.130526192220051,
    0.793353340291235, -0.60876142900872,
    0.608761429008721, -0.793353340291235,
    0.130526192220052, -0.99144486137381,
    -0.130526192220052, -0.99144486137381,
    -0.608761429008721, -0.793353340291235,
    -0.793353340291235, -0.608761429008721,
    -0.99144486137381,  -0.130526192220052,
    -0.99144486137381,   0.130526192220051,
    -0.793353340291235,  0.608761429008721,
    -0.608761429008721,  0.793353340291235,
    -0.130526192220052,  0.99144486137381,
];

#[rustfmt::skip]
const GRAD3_SRC: &[f64] = &[
    2.22474487139,       2.22474487139,      -1.0,                 0.0,
    2.22474487139,       2.22474487139,       1.0,                 0.0,
    3.0862664687972017,  1.1721513422464978,  0.0,                 0.0,
    1.1721513422464978,  3.0862664687972017,  0.0,                 0.0,
    -2.22474487139,       2.22474487139,      -1.0,                 0.0,
    -2.22474487139,       2.22474487139,       1.0,                 0.0,
    -1.1721513422464978,  3.0862664687972017,  0.0,                 0.0,
    -3.0862664687972017,  1.1721513422464978,  0.0,                 0.0,
    -1.0,                -2.22474487139,      -2.22474487139,       0.0,
    1.0,                -2.22474487139,      -2.22474487139,       0.0,
    0.0,                -3.0862664687972017, -1.1721513422464978,  0.0,
    0.0,                -1.1721513422464978, -3.0862664687972017,  0.0,
    -1.0,                -2.22474487139,       2.22474487139,       0.0,
    1.0,                -2.22474487139,       2.22474487139,       0.0,
    0.0,                -1.1721513422464978,  3.0862664687972017,  0.0,
    0.0,                -3.0862664687972017,  1.1721513422464978,  0.0,
    //--------------------------------------------------------------------//
    -2.22474487139,      -2.22474487139,      -1.0,                 0.0,
    -2.22474487139,      -2.22474487139,       1.0,                 0.0,
    -3.0862664687972017, -1.1721513422464978,  0.0,                 0.0,
    -1.1721513422464978, -3.0862664687972017,  0.0,                 0.0,
    -2.22474487139,      -1.0,                -2.22474487139,       0.0,
    -2.22474487139,       1.0,                -2.22474487139,       0.0,
    -1.1721513422464978,  0.0,                -3.0862664687972017,  0.0,
    -3.0862664687972017,  0.0,                -1.1721513422464978,  0.0,
    -2.22474487139,      -1.0,                 2.22474487139,       0.0,
    -2.22474487139,       1.0,                 2.22474487139,       0.0,
    -3.0862664687972017,  0.0,                 1.1721513422464978,  0.0,
    -1.1721513422464978,  0.0,                 3.0862664687972017,  0.0,
    -1.0,                 2.22474487139,      -2.22474487139,       0.0,
    1.0,                 2.22474487139,      -2.22474487139,       0.0,
    0.0,                 1.1721513422464978, -3.0862664687972017,  0.0,
    0.0,                 3.0862664687972017, -1.1721513422464978,  0.0,
    -1.0,                 2.22474487139,       2.22474487139,       0.0,
    1.0,                 2.22474487139,       2.22474487139,       0.0,
    0.0,                 3.0862664687972017,  1.1721513422464978,  0.0,
    0.0,                 1.1721513422464978,  3.0862664687972017,  0.0,
    2.22474487139,      -2.22474487139,      -1.0,                 0.0,
    2.22474487139,      -2.22474487139,       1.0,                 0.0,
    1.1721513422464978, -3.0862664687972017,  0.0,                 0.0,
    3.0862664687972017, -1.1721513422464978,  0.0,                 0.0,
    2.22474487139,      -1.0,                -2.22474487139,       0.0,
    2.22474487139,       1.0,                -2.22474487139,       0.0,
    3.0862664687972017,  0.0,                -1.1721513422464978,  0.0,
    1.1721513422464978,  0.0,                -3.0862664687972017,  0.0,
    2.22474487139,      -1.0,                 2.22474487139,       0.0,
    2.22474487139,       1.0,                 2.22474487139,       0.0,
    1.1721513422464978,  0.0,                 3.0862664687972017,  0.0,
    3.0862664687972017,  0.0,                 1.1721513422464978,  0.0,
];

#[rustfmt::skip]
const GRAD4_SRC: &[f64] = &[
    -0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537,    0.5794684678643381,
    -0.7504883828755602,   -0.4004672082940195,    0.15296486218853164,   0.5029860367700724,
    -0.7504883828755602,    0.15296486218853164,  -0.4004672082940195,    0.5029860367700724,
    -0.8828161875373585,    0.08164729285680945,   0.08164729285680945,   0.4553054119602712,
    -0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945,   0.8828161875373585,
    -0.5029860367700724,   -0.15296486218853164,   0.4004672082940195,    0.7504883828755602,
    -0.5029860367700724,    0.4004672082940195,   -0.15296486218853164,   0.7504883828755602,
    -0.5794684678643381,    0.3239847771997537,    0.3239847771997537,    0.6740059517812944,
    -0.6740059517812944,   -0.3239847771997537,    0.5794684678643381,   -0.3239847771997537,
    -0.7504883828755602,   -0.4004672082940195,    0.5029860367700724,    0.15296486218853164,
    -0.7504883828755602,    0.15296486218853164,   0.5029860367700724,   -0.4004672082940195,
    -0.8828161875373585,    0.08164729285680945,   0.4553054119602712,    0.08164729285680945,
    -0.4553054119602712,   -0.08164729285680945,   0.8828161875373585,   -0.08164729285680945,
    -0.5029860367700724,   -0.15296486218853164,   0.7504883828755602,    0.4004672082940195,
    -0.5029860367700724,    0.4004672082940195,    0.7504883828755602,   -0.15296486218853164,
    -0.5794684678643381,    0.3239847771997537,    0.6740059517812944,    0.3239847771997537,
    -0.6740059517812944,    0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537,
    -0.7504883828755602,    0.5029860367700724,   -0.4004672082940195,    0.15296486218853164,
    -0.7504883828755602,    0.5029860367700724,    0.15296486218853164,  -0.4004672082940195,
    -0.8828161875373585,    0.4553054119602712,    0.08164729285680945,   0.08164729285680945,
    -0.4553054119602712,    0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945,
    -0.5029860367700724,    0.7504883828755602,   -0.15296486218853164,   0.4004672082940195,
    -0.5029860367700724,    0.7504883828755602,    0.4004672082940195,   -0.15296486218853164,
    -0.5794684678643381,    0.6740059517812944,    0.3239847771997537,    0.3239847771997537,
    0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537,
    0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195,    0.15296486218853164,
    0.5029860367700724,   -0.7504883828755602,    0.15296486218853164,  -0.4004672082940195,
    0.4553054119602712,   -0.8828161875373585,    0.08164729285680945,   0.08164729285680945,
    0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945,
    0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164,   0.4004672082940195,
    0.7504883828755602,   -0.5029860367700724,    0.4004672082940195,   -0.15296486218853164,
    0.6740059517812944,   -0.5794684678643381,    0.3239847771997537,    0.3239847771997537,
    //------------------------------------------------------------------------------------------//
    -0.753341017856078,    -0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624,
    -0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301,    0.12128480194602098,
    -0.7821684431180708,   -0.4321472685365301,    0.12128480194602098,  -0.4321472685365301,
    -0.7821684431180708,    0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301,
    -0.8586508742123365,   -0.508629699630796,     0.044802370851755174,  0.044802370851755174,
    -0.8586508742123365,    0.044802370851755174, -0.508629699630796,     0.044802370851755174,
    -0.8586508742123365,    0.044802370851755174,  0.044802370851755174, -0.508629699630796,
    -0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842,
    -0.37968289875261624,  -0.753341017856078,    -0.37968289875261624,  -0.37968289875261624,
    -0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301,    0.12128480194602098,
    -0.4321472685365301,   -0.7821684431180708,    0.12128480194602098,  -0.4321472685365301,
    0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301,
    -0.508629699630796,    -0.8586508742123365,    0.044802370851755174,  0.044802370851755174,
    0.044802370851755174, -0.8586508742123365,   -0.508629699630796,     0.044802370851755174,
    0.044802370851755174, -0.8586508742123365,    0.044802370851755174, -0.508629699630796,
    -0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842,
    -0.37968289875261624,  -0.37968289875261624,  -0.753341017856078,    -0.37968289875261624,
    -0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708,    0.12128480194602098,
    -0.4321472685365301,    0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301,
    0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301,
    -0.508629699630796,     0.044802370851755174, -0.8586508742123365,    0.044802370851755174,
    0.044802370851755174, -0.508629699630796,    -0.8586508742123365,    0.044802370851755174,
    0.044802370851755174,  0.044802370851755174, -0.8586508742123365,   -0.508629699630796,
    -0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842,
    -0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624,  -0.753341017856078,
    -0.4321472685365301,   -0.4321472685365301,    0.12128480194602098,  -0.7821684431180708,
    -0.4321472685365301,    0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708,
    0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708,
    -0.508629699630796,     0.044802370851755174,  0.044802370851755174, -0.8586508742123365,
    0.044802370851755174, -0.508629699630796,     0.044802370851755174, -0.8586508742123365,
    0.044802370851755174,  0.044802370851755174, -0.508629699630796,    -0.8586508742123365,
    -0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062,
    -0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537,    0.5794684678643381,
    -0.4004672082940195,   -0.7504883828755602,    0.15296486218853164,   0.5029860367700724,
    0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195,    0.5029860367700724,
    0.08164729285680945,  -0.8828161875373585,    0.08164729285680945,   0.4553054119602712,
    -0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945,   0.8828161875373585,
    -0.15296486218853164,  -0.5029860367700724,    0.4004672082940195,    0.7504883828755602,
    0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164,   0.7504883828755602,
    0.3239847771997537,   -0.5794684678643381,    0.3239847771997537,    0.6740059517812944,
    -0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944,    0.5794684678643381,
    -0.4004672082940195,    0.15296486218853164,  -0.7504883828755602,    0.5029860367700724,
    0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602,    0.5029860367700724,
    0.08164729285680945,   0.08164729285680945,  -0.8828161875373585,    0.4553054119602712,
    -0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712,    0.8828161875373585,
    -0.15296486218853164,   0.4004672082940195,   -0.5029860367700724,    0.7504883828755602,
    0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724,    0.7504883828755602,
    0.3239847771997537,    0.3239847771997537,   -0.5794684678643381,    0.6740059517812944,
    -0.3239847771997537,   -0.6740059517812944,    0.5794684678643381,   -0.3239847771997537,
    -0.4004672082940195,   -0.7504883828755602,    0.5029860367700724,    0.15296486218853164,
    0.15296486218853164,  -0.7504883828755602,    0.5029860367700724,   -0.4004672082940195,
    0.08164729285680945,  -0.8828161875373585,    0.4553054119602712,    0.08164729285680945,
    -0.08164729285680945,  -0.4553054119602712,    0.8828161875373585,   -0.08164729285680945,
    -0.15296486218853164,  -0.5029860367700724,    0.7504883828755602,    0.4004672082940195,
    0.4004672082940195,   -0.5029860367700724,    0.7504883828755602,   -0.15296486218853164,
    0.3239847771997537,   -0.5794684678643381,    0.6740059517812944,    0.3239847771997537,
    -0.3239847771997537,   -0.3239847771997537,    0.5794684678643381,   -0.6740059517812944,
    -0.4004672082940195,    0.15296486218853164,   0.5029860367700724,   -0.7504883828755602,
    0.15296486218853164,  -0.4004672082940195,    0.5029860367700724,   -0.7504883828755602,
    0.08164729285680945,   0.08164729285680945,   0.4553054119602712,   -0.8828161875373585,
    -0.08164729285680945,  -0.08164729285680945,   0.8828161875373585,   -0.4553054119602712,
    -0.15296486218853164,   0.4004672082940195,    0.7504883828755602,   -0.5029860367700724,
    0.4004672082940195,   -0.15296486218853164,   0.7504883828755602,   -0.5029860367700724,
    0.3239847771997537,    0.3239847771997537,    0.6740059517812944,   -0.5794684678643381,
    -0.3239847771997537,    0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537,
    -0.4004672082940195,    0.5029860367700724,   -0.7504883828755602,    0.15296486218853164,
    0.15296486218853164,   0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195,
    0.08164729285680945,   0.4553054119602712,   -0.8828161875373585,    0.08164729285680945,
    -0.08164729285680945,   0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945,
    -0.15296486218853164,   0.7504883828755602,   -0.5029860367700724,    0.4004672082940195,
    0.4004672082940195,    0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164,
    0.3239847771997537,    0.6740059517812944,   -0.5794684678643381,    0.3239847771997537,
    -0.3239847771997537,    0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944,
    -0.4004672082940195,    0.5029860367700724,    0.15296486218853164,  -0.7504883828755602,
    0.15296486218853164,   0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602,
    0.08164729285680945,   0.4553054119602712,    0.08164729285680945,  -0.8828161875373585,
    -0.08164729285680945,   0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712,
    -0.15296486218853164,   0.7504883828755602,    0.4004672082940195,   -0.5029860367700724,
    0.4004672082940195,    0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724,
    0.3239847771997537,    0.6740059517812944,    0.3239847771997537,   -0.5794684678643381,
    0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537,
    0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602,    0.15296486218853164,
    0.5029860367700724,    0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195,
    0.4553054119602712,    0.08164729285680945,  -0.8828161875373585,    0.08164729285680945,
    0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945,
    0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724,    0.4004672082940195,
    0.7504883828755602,    0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164,
    0.6740059517812944,    0.3239847771997537,   -0.5794684678643381,    0.3239847771997537,
    0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944,
    0.5029860367700724,   -0.4004672082940195,    0.15296486218853164,  -0.7504883828755602,
    0.5029860367700724,    0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602,
    0.4553054119602712,    0.08164729285680945,   0.08164729285680945,  -0.8828161875373585,
    0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712,
    0.7504883828755602,   -0.15296486218853164,   0.4004672082940195,   -0.5029860367700724,
    0.7504883828755602,    0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724,
    0.6740059517812944,    0.3239847771997537,    0.3239847771997537,   -0.5794684678643381,
    0.03381941603233842,   0.03381941603233842,   0.03381941603233842,   0.9982828964265062,
    -0.044802370851755174, -0.044802370851755174,  0.508629699630796,     0.8586508742123365,
    -0.044802370851755174,  0.508629699630796,    -0.044802370851755174,  0.8586508742123365,
    -0.12128480194602098,   0.4321472685365301,    0.4321472685365301,    0.7821684431180708,
    0.508629699630796,    -0.044802370851755174, -0.044802370851755174,  0.8586508742123365,
    0.4321472685365301,   -0.12128480194602098,   0.4321472685365301,    0.7821684431180708,
    0.4321472685365301,    0.4321472685365301,   -0.12128480194602098,   0.7821684431180708,
    0.37968289875261624,   0.37968289875261624,   0.37968289875261624,   0.753341017856078,
    0.03381941603233842,   0.03381941603233842,   0.9982828964265062,    0.03381941603233842,
    -0.044802370851755174,  0.044802370851755174,  0.8586508742123365,    0.508629699630796,
    -0.044802370851755174,  0.508629699630796,     0.8586508742123365,   -0.044802370851755174,
    -0.12128480194602098,   0.4321472685365301,    0.7821684431180708,    0.4321472685365301,
    0.508629699630796,    -0.044802370851755174,  0.8586508742123365,   -0.044802370851755174,
    0.4321472685365301,   -0.12128480194602098,   0.7821684431180708,    0.4321472685365301,
    0.4321472685365301,    0.4321472685365301,    0.7821684431180708,   -0.12128480194602098,
    0.37968289875261624,   0.37968289875261624,   0.753341017856078,     0.37968289875261624,
    0.03381941603233842,   0.9982828964265062,    0.03381941603233842,   0.03381941603233842,
    -0.044802370851755174,  0.8586508742123365,   -0.044802370851755174,  0.508629699630796,
    -0.044802370851755174,  0.8586508742123365,    0.508629699630796,    -0.044802370851755174,
    -0.12128480194602098,   0.7821684431180708,    0.4321472685365301,    0.4321472685365301,
    0.508629699630796,     0.8586508742123365,   -0.044802370851755174, -0.044802370851755174,
    0.4321472685365301,    0.7821684431180708,   -0.12128480194602098,   0.4321472685365301,
    0.4321472685365301,    0.7821684431180708,    0.4321472685365301,   -0.12128480194602098,
    0.37968289875261624,   0.753341017856078,     0.37968289875261624,   0.37968289875261624,
    0.9982828964265062,    0.03381941603233842,   0.03381941603233842,   0.03381941603233842,
    0.8586508742123365,   -0.044802370851755174, -0.044802370851755174,  0.508629699630796,
    0.8586508742123365,   -0.044802370851755174,  0.508629699630796,    -0.044802370851755174,
    0.7821684431180708,   -0.12128480194602098,   0.4321472685365301,    0.4321472685365301,
    0.8586508742123365,    0.508629699630796,    -0.044802370851755174, -0.044802370851755174,
    0.7821684431180708,    0.4321472685365301,   -0.12128480194602098,   0.4321472685365301,
    0.7821684431180708,    0.4321472685365301,    0.4321472685365301,   -0.12128480194602098,
    0.753341017856078,     0.37968289875261624,   0.37968289875261624,   0.37968289875261624,
];

#[rustfmt::skip]
const LOOKUP_4D_VERTEX_CODES: &[&[u8]] = &[
    &[0x15, 0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA],
    &[0x01, 0x05, 0x11, 0x15, 0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA],
    &[0x01, 0x15, 0x16, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA],
    &[0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA],
    &[0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA],
    &[0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA, 0xAB],
    &[0x04, 0x05, 0x14, 0x15, 0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA],
    &[0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA],
    &[0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA],
    &[0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x6A, 0x9A, 0xAA, 0xAB],
    &[0x04, 0x15, 0x19, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA, 0xAE],
    &[0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x5E, 0x6A, 0x9A, 0xAA, 0xAE],
    &[0x05, 0x15, 0x1A, 0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x5E, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x11, 0x15, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA],
    &[0x11, 0x15, 0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA],
    &[0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA, 0xAB],
    &[0x14, 0x15, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA6, 0xA9, 0xAA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x15, 0x16, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x6B, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x14, 0x15, 0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xAA],
    &[0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x6B, 0x9A, 0xAA, 0xAB],
    &[0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x19, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x9A, 0xAA, 0xAE],
    &[0x15, 0x1A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x10, 0x11, 0x14, 0x15, 0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA],
    &[0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA],
    &[0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x67, 0x6A, 0xA6, 0xAA, 0xAB],
    &[0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA],
    &[0x15, 0x16, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x6B, 0xA6, 0xAA, 0xAB],
    &[0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA],
    &[0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0xAA, 0xAB],
    &[0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x6D, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x19, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x6E, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0xAA, 0xAE],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x10, 0x15, 0x25, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA, 0xBA],
    &[0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0x76, 0xA6, 0xAA, 0xBA],
    &[0x11, 0x15, 0x26, 0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x67, 0x6A, 0x76, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB],
    &[0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x25, 0x55, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x25, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA6, 0xAA, 0xBA],
    &[0x15, 0x26, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB],
    &[0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0x79, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x25, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x25, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xAA, 0xBA],
    &[0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB],
    &[0x14, 0x15, 0x29, 0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0x6D, 0x79, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE],
    &[0x15, 0x29, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE],
    &[0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x7A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF],
    &[0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x41, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA],
    &[0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA],
    &[0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x44, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x45, 0x46, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB],
    &[0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA],
    &[0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x9A, 0x9B, 0xAA, 0xAB],
    &[0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x45, 0x49, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE],
    &[0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x99, 0x9A, 0x9E, 0xAA, 0xAE],
    &[0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x50, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x51, 0x55, 0x56, 0x59, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB],
    &[0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB],
    &[0x54, 0x55, 0x56, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE],
    &[0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE],
    &[0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA],
    &[0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xA7, 0xAA, 0xAB],
    &[0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA],
    &[0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA],
    &[0x15, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAD, 0xAE],
    &[0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x51, 0x55, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA],
    &[0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA],
    &[0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x6A, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB],
    &[0x54, 0x55, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA],
    &[0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB],
    &[0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA],
    &[0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB],
    &[0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE],
    &[0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE],
    &[0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE],
    &[0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF],
    &[0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA],
    &[0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA],
    &[0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA],
    &[0x45, 0x46, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB],
    &[0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA],
    &[0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xAA, 0xAB],
    &[0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x95, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE],
    &[0x45, 0x49, 0x55, 0x59, 0x5A, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE],
    &[0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xAA, 0xAE],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA],
    &[0x51, 0x52, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xEA],
    &[0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA],
    &[0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xEA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB],
    &[0x54, 0x55, 0x58, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE],
    &[0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA],
    &[0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA],
    &[0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB],
    &[0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE],
    &[0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE],
    &[0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAE],
    &[0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xAF],
    &[0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA],
    &[0x51, 0x55, 0x61, 0x65, 0x66, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA],
    &[0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB],
    &[0x54, 0x55, 0x64, 0x65, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA],
    &[0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB],
    &[0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA],
    &[0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA],
    &[0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xBB],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE],
    &[0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE],
    &[0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA],
    &[0x40, 0x45, 0x51, 0x54, 0x55, 0x85, 0x91, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA, 0xEA],
    &[0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xD6, 0xEA],
    &[0x41, 0x45, 0x51, 0x55, 0x56, 0x86, 0x92, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB, 0xD6, 0xEA, 0xEB],
    &[0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xEA],
    &[0x45, 0x55, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xDA, 0xEA],
    &[0x45, 0x55, 0x56, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xDA, 0xEA],
    &[0x45, 0x55, 0x56, 0x86, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB],
    &[0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xD9, 0xEA],
    &[0x45, 0x55, 0x59, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xDA, 0xEA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xDA, 0xEA],
    &[0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB],
    &[0x44, 0x45, 0x54, 0x55, 0x59, 0x89, 0x95, 0x98, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE, 0xD9, 0xEA, 0xEE],
    &[0x45, 0x55, 0x59, 0x89, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE],
    &[0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xDA, 0xEA, 0xEF],
    &[0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x51, 0x55, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA],
    &[0x51, 0x55, 0x56, 0x91, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA],
    &[0x51, 0x55, 0x56, 0x92, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB],
    &[0x54, 0x55, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA],
    &[0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB],
    &[0x54, 0x55, 0x59, 0x94, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA],
    &[0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB],
    &[0x54, 0x55, 0x59, 0x95, 0x98, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE],
    &[0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE],
    &[0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE],
    &[0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xEA, 0xEF],
    &[0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xE5, 0xEA],
    &[0x51, 0x55, 0x65, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA],
    &[0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB],
    &[0x54, 0x55, 0x65, 0x94, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA],
    &[0x55, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x94, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x55, 0x56, 0x59, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA],
    &[0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xEA, 0xEB],
    &[0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE],
    &[0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE],
    &[0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE],
    &[0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xEA],
    &[0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA1, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA, 0xE5, 0xEA, 0xFA],
    &[0x51, 0x55, 0x65, 0x95, 0xA1, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA],
    &[0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA],
    &[0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xE6, 0xEA, 0xFB],
    &[0x54, 0x55, 0x65, 0x95, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA],
    &[0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA],
    &[0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA],
    &[0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xEA, 0xFB],
    &[0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA],
    &[0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA],
    &[0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA],
    &[0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xEA],
    &[0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xE9, 0xEA, 0xFE],
    &[0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA, 0xFE],
    &[0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA],
    &[0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA, 0xEA],
];
