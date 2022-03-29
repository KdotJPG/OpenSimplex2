#ifndef __bcc_noise_4_hlsl_
#define __bcc_noise_4_hlsl_

/////////////// K.jpg's Simplex-like Re-oriented 4-Point BCC Noise ///////////////
//////////////////// Output: float4(dF/dx, dF/dy, dF/dz, value) ////////////////////

// Inspired by Stefan Gustavson's noise
float4 permute(float4 t) {
	return t * (t * 34.0 + 133.0);
}

// https://stackoverflow.com/questions/7610631/glsl-mod-vs-hlsl-fmod
float mod(float x, float y)
{
	return x - y * floor(x / y);
}

float2 mod(float2 x, float2 y)
{
	return x - y * floor(x / y);
}

float3 mod(float3 x, float3 y)
{
	return x - y * floor(x / y);
}

float4 mod(float4 x, float4 y)
{
	return x - y * floor(x / y);
}

// Gradient set is a normalized expanded rhombic dodecahedron
float3 grad(float hash) {
	// Random vertex of a cube, +/- 1 each
	float3 cube = mod(floor(hash / float3(1.0, 2.0, 4.0)), 2.0) * 2.0 - 1.0;

	// Random edge of the three edges connected to that vertex
	// Also a cuboctahedral vertex
	// And corresponds to the face of its dual, the rhombic dodecahedron
	float3 cuboct = cube;

	int index = int(hash / 16.0);

	if (index == 0)
		cuboct.x = 0.0;
	else if (index == 1)
		cuboct.y = 0.0;
	else
		cuboct.z = 0.0;

	// In a funky way, pick one of the four points on the rhombic face
	float type = mod(floor(hash / 8.0), 2.0);
	float3 rhomb = (1.0 - type) * cube + type * (cuboct + cross(cube, cuboct));

	// Expand it so that the new edges are the same length
	// as the existing ones
	float3 grad = cuboct * 1.22474487139 + rhomb;

	// To make all gradients the same length, we only need to shorten the
	// second type of vector. We also put in the whole noise scale constant.
	// The compiler should reduce it into the existing floats. I think.
	grad *= (1.0 - 0.042942436724648037 * type) * 32.80201376986577;

	return grad;
}

// BCC lattice split up into 2 cube lattices
float4 openSimplex2Base(float3 X) {
	// First half-lattice, closest edge
	float3 v1 = round(X);
	float3 d1 = X - v1;
	float3 score1 = abs(d1);
	float3 dir1 = step(max(score1.yzx, score1.zxy), score1);
	float3 v2 = v1 + dir1 * sign(d1);
	float3 d2 = X - v2;

	// Second half-lattice, closest edge
	float3 X2 = X + 144.5;
	float3 v3 = round(X2);
	float3 d3 = X2 - v3;
	float3 score2 = abs(d3);
	float3 dir2 = step(max(score2.yzx, score2.zxy), score2);
	float3 v4 = v3 + dir2 * sign(d3);
	float3 d4 = X2 - v4;

	// Gradient hashes for the four points, two from each half-lattice
	float4 hashes = permute(mod(float4(v1.x, v2.x, v3.x, v4.x), 289.0));
	hashes = permute(mod(hashes + float4(v1.y, v2.y, v3.y, v4.y), 289.0));
	hashes = mod(permute(mod(hashes + float4(v1.z, v2.z, v3.z, v4.z), 289.0)), 48.0);

	// Gradient extrapolations & kernel function
	float4 a = max(0.5 - float4(dot(d1, d1), dot(d2, d2), dot(d3, d3), dot(d4, d4)), 0.0);
	float4 aa = a * a; float4 aaaa = aa * aa;
	float3 g1 = grad(hashes.x); float3 g2 = grad(hashes.y);
	float3 g3 = grad(hashes.z); float3 g4 = grad(hashes.w);
	float4 extrapolations = float4(dot(d1, g1), dot(d2, g2), dot(d3, g3), dot(d4, g4));

	float4x3 derivativeMatrix = { d1, d2, d3, d4 };
	float4x3 gradientMatrix = { g1, g2, g3, g4 };

	// Derivatives of the noise
	float3 derivative = -8.0 * mul(aa * a * extrapolations, derivativeMatrix)
		+ mul(aaaa, gradientMatrix);

	// Return it all as a float4
	return float4(derivative, dot(aaaa, extrapolations));
}

// Use this if you don't want Z to look different from X and Y
float4 openSimplex2_Conventional(float3 X) {
	// Rotate around the main diagonal. Not a skew transform.
	float4 result = openSimplex2Base(dot(X, float3(2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)) - X);
	return float4(dot(result.xyz, float3(2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)) - result.xyz, result.w);
}

// Use this if you want to show X and Y in a plane, then use Z for time, vertical, etc.
float4 openSimplex2_ImproveXY(float3 X) {
	// Rotate so Z points down the main diagonal. Not a skew transform.
	float3x3 orthonormalMap = {
		0.788675134594813, -0.211324865405187, -0.577350269189626,
			-0.211324865405187, 0.788675134594813, -0.577350269189626,
			0.577350269189626, 0.577350269189626, 0.577350269189626 };

	float4 result = openSimplex2Base(mul(X, orthonormalMap));
	return float4(mul(orthonormalMap, result.xyz), result.w);
}

//////////////////////////////// End noise code ////////////////////////////////

#endif
