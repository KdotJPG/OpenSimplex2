use std::ffi::{c_double, c_float, c_longlong};

use crate::{fast, smooth};

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise2(seed: c_longlong, x: c_double, y: c_double) -> c_float {
    fast::noise2(seed, x, y)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise2_ImproveX(
    seed: c_longlong,
    x: c_double,
    y: c_double,
) -> c_float {
    fast::noise2_ImproveX(seed, x, y)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise3_ImproveXY(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
) -> c_float {
    fast::noise3_ImproveXY(seed, x, y, z)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise3_ImproveXZ(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
) -> c_float {
    fast::noise3_ImproveXZ(seed, x, y, z)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise3_Fallback(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
) -> c_float {
    fast::noise3_Fallback(seed, x, y, z)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise4_ImproveXYZ_ImproveXY(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    fast::noise4_ImproveXYZ_ImproveXY(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise4_ImproveXYZ_ImproveXZ(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    fast::noise4_ImproveXYZ_ImproveXZ(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise4_ImproveXYZ(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    fast::noise4_ImproveXYZ(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise4_ImproveXY_ImproveZW(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    fast::noise4_ImproveXY_ImproveZW(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_fast_noise4_Fallback(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    fast::noise4_Fallback(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise2(
    seed: c_longlong,
    x: c_double,
    y: c_double,
) -> c_float {
    smooth::noise2(seed, x, y)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise2_ImproveX(
    seed: c_longlong,
    x: c_double,
    y: c_double,
) -> c_float {
    smooth::noise2_ImproveX(seed, x, y)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise3_ImproveXY(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
) -> c_float {
    smooth::noise3_ImproveXY(seed, x, y, z)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise3_ImproveXZ(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
) -> c_float {
    smooth::noise3_ImproveXZ(seed, x, y, z)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise3_Fallback(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
) -> c_float {
    smooth::noise3_Fallback(seed, x, y, z)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise4_ImproveXYZ_ImproveXY(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    smooth::noise4_ImproveXYZ_ImproveXY(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise4_ImproveXYZ_ImproveXZ(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    smooth::noise4_ImproveXYZ_ImproveXZ(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise4_ImproveXYZ(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    smooth::noise4_ImproveXYZ(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise4_ImproveXY_ImproveZW(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    smooth::noise4_ImproveXY_ImproveZW(seed, x, y, z, w)
}

#[no_mangle]
pub extern "C" fn opensimplex2_smooth_noise4_Fallback(
    seed: c_longlong,
    x: c_double,
    y: c_double,
    z: c_double,
    w: c_double,
) -> c_float {
    smooth::noise4_Fallback(seed, x, y, z, w)
}
