#pragma once
#ifndef OPENSIMPLEX2_H

extern float opensimplex2_fast_noise2(long long seed, double x, double y);
extern float opensimplex2_fast_noise2_ImproveX(long long seed, double x, double y);
extern float opensimplex2_fast_noise3_ImproveXY(long long seed, double x, double y, double z);
extern float opensimplex2_fast_noise3_ImproveXZ(long long seed, double x, double y, double z);
extern float opensimplex2_fast_noise3_Fallback(long long seed, double x, double y, double z);
extern float opensimplex2_fast_noise4_ImproveXYZ_ImproveXY(long long seed, double x, double y, double z, double w);
extern float opensimplex2_fast_noise4_ImproveXYZ_ImproveXZ(long long seed, double x, double y, double z, double w);
extern float opensimplex2_fast_noise4_ImproveXYZ(long long seed, double x, double y, double z, double w);
extern float opensimplex2_fast_noise4_ImproveXY_ImproveZW(long long seed, double x, double y, double z, double w);
extern float opensimplex2_fast_noise4_Fallback(long long seed, double x, double y, double z, double w);
extern float opensimplex2_smooth_noise2(long long seed, double x, double y);
extern float opensimplex2_smooth_noise2_ImproveX(long long seed, double x, double y);
extern float opensimplex2_smooth_noise3_ImproveXY(long long seed, double x, double y, double z);
extern float opensimplex2_smooth_noise3_ImproveXZ(long long seed, double x, double y, double z);
extern float opensimplex2_smooth_noise3_Fallback(long long seed, double x, double y, double z);
extern float opensimplex2_smooth_noise4_ImproveXYZ_ImproveXY(long long seed, double x, double y, double z, double w);
extern float opensimplex2_smooth_noise4_ImproveXYZ_ImproveXZ(long long seed, double x, double y, double z, double w);
extern float opensimplex2_smooth_noise4_ImproveXYZ(long long seed, double x, double y, double z, double w);
extern float opensimplex2_smooth_noise4_ImproveXY_ImproveZW(long long seed, double x, double y, double z, double w);
extern float opensimplex2_smooth_noise4_Fallback(long long seed, double x, double y, double z, double w);

#endif
