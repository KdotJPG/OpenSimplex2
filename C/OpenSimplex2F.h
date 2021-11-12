#ifndef OpenSimplex2F_h__
#define OpenSimplex2F_h__

/**
 * K.jpg's OpenSimplex 2, faster variant
 *
 * - 2D is standard simplex implemented using a lookup table.
 * - 3D is "Re-oriented 4-point BCC noise" which constructs a
 *   congruent BCC lattice in a much different way than usual.
 * - 4D constructs the lattice as a union of five copies of its
 *   reciprocal. It successively finds the closest point on each.
 *
 * Multiple versions of each function are provided. See the
 * documentation above each, for more info.
 *
 * Ported from Java to C by Stephen M. Cameron
 *
 */

#if ((__GNUC_STDC_INLINE__) || (__STDC_VERSION__ >= 199901L))
	#include <stdint.h>
	#define INLINE inline
#elif (defined (_MSC_VER) || defined (__GNUC_GNU_INLINE__))
	#include <stdint.h>
	#define INLINE __inline
#else
	/* ANSI C doesn't have inline or stdint.h. */
	#define INLINE
#endif

#ifdef __cplusplus
	extern "C" {
#endif

struct OpenSimplex2F_context;

/* Allocate and initialize OpenSimplex2F context */
int OpenSimplex2F(int64_t seed, struct OpenSimplex2F_context **ctx);

/* Free OpenSimplex2F context */
void OpenSimplex2F_free(struct OpenSimplex2F_context * ctx);

/* Free singleton lattice point data */
void OpenSimplex2F_shutdown(void);

/**
 * 2D Simplex noise, standard lattice orientation.
 */
double OpenSimplex2F_noise2(struct OpenSimplex2F_context *ctx, double x, double y);

/**
 * 2D Simplex noise, with Y pointing down the main diagonal.
 * Might be better for a 2D sandbox style game, where Y is vertical.
 * Probably slightly less optimal for heightmaps or continent maps.
 */
double OpenSimplex2F_noise2_XBeforeY(struct OpenSimplex2F_context *ctx, double x, double y);

/**
 * 3D Re-oriented 4-point BCC noise, classic orientation.
 * Proper substitute for 3D Simplex in light of Forbidden Formulae.
 * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
 */
double OpenSimplex2F_noise3_Classic(struct OpenSimplex2F_context *ctx, double x, double y, double z);

/**
 * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Y).
 * Recommended for 3D terrain and time-varied animations.
 * The Z coordinate should always be the "different" coordinate in your use case.
 * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
 * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
 * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
 */
double OpenSimplex2F_noise3_XYBeforeZ(struct OpenSimplex2F_context *ctx, double x, double y, double z);

/**
 * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Z).
 * Recommended for 3D terrain and time-varied animations.
 * The Y coordinate should always be the "different" coordinate in your use case.
 * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
 * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
 * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
 */
double OpenSimplex2F_noise3_XZBeforeY(struct OpenSimplex2F_context *ctx, double x, double y, double z);

/**
 * 4D OpenSimplex2F noise, classic lattice orientation.
 */
double OpenSimplex2F_noise4_Classic(struct OpenSimplex2F_context *ctx, double x, double y, double z, double w);

/**
 * 4D OpenSimplex2F noise, with XY and ZW forming orthogonal triangular-based planes.
 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
 * Recommended for noise(x, y, sin(time), cos(time)) trick.
 */
double OpenSimplex2F_noise4_XYBeforeZW(struct OpenSimplex2F_context *ctx, double x, double y, double z, double w);

/**
 * 4D OpenSimplex2F noise, with XZ and YW forming orthogonal triangular-based planes.
 * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
 */
double OpenSimplex2F_noise4_XZBeforeYW(struct OpenSimplex2F_context *ctx, double x, double y, double z, double w);

/**
 * 4D OpenSimplex2F noise, with XYZ oriented like noise3_Classic,
 * and W for an extra degree of freedom. W repeats eventually.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 */
double OpenSimplex2F_noise4_XYZBeforeW(struct OpenSimplex2F_context *ctx, double x, double y, double z, double w);

#ifdef __cplusplus
	}
#endif

#endif

