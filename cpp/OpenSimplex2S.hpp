/**
 * Ported from https://github.com/KdotJPG/OpenSimplex2/blob/master/java/OpenSimplex2S.java
 * Probably not best implementation of static initialization.
 * Also changed some code to use fixed c-style arrays, to avoid using of std library.
 */

/**
 * K.jpg's OpenSimplex 2, smooth variant ("SuperSimplex")
 *
 * - 2D is standard simplex, modified to support larger kernels.
 *   Implemented using a lookup table.
 * - 3D is "Re-oriented 8-point BCC noise" which constructs a
 *   congruent BCC lattice in a much different way than usual.
 * - 4D uses a naïve pregenerated lookup table, and averages out
 *   to the expected performance.
 *
 * Multiple versions of each function are provided. See the
 * documentation above each, for more info.
 */
class OpenSimplex2S
{
  struct Grad2
  {
    double dx {};
    double dy {};
  };

  struct Grad3
  {
    double dx {};
    double dy {};
    double dz {};
  };

  struct Grad4
  {
    double dx {};
    double dy {};
    double dz {};
    double dw {};
  };

  struct LatticePoint2D
  {
    int xsv {};
    int ysv {};
    double dx {};
    double dy {};

    LatticePoint2D();
    LatticePoint2D(int xsv, int ysv);
  };

  struct LatticePoint3D
  {
    double dxr {};
    double dyr {};
    double dzr {};
    int xrv {};
    int yrv {};
    int zrv {};
    LatticePoint3D* nextOnFailure {};
    LatticePoint3D* nextOnSuccess {};
    LatticePoint3D();
    LatticePoint3D(int xrv, int yrv, int zrv, int lattice);
  };

  struct LatticePoint4D
  {
    int xsv {};
    int ysv {};
    int zsv {};
    int wsv {};
    double dx {};
    double dy {};
    double dz {};
    double dw {};

    LatticePoint4D();
    LatticePoint4D(int xsv, int ysv, int zsv, int wsv);
  };

  static const int PSIZE = 2048;
  static const int PMASK = 2047;

  constexpr static const double N2 = 0.05481866495625118;
  constexpr static const double N3 = 0.2781926117527186;
  constexpr static const double N4 = 0.11127401889945551;

  static Grad2 GRADIENTS_2D[PSIZE];
  static Grad3 GRADIENTS_3D[PSIZE];
  static Grad4 GRADIENTS_4D[PSIZE];

  static LatticePoint2D LOOKUP_2D[8 * 4];
  static LatticePoint3D LOOKUP_3D[8];
  static LatticePoint4D LOOKUP_4D[256][20];
  static unsigned char LOOKUP_4D_SIZE[256];

  short perm[PSIZE];
  Grad2 permGrad2[PSIZE];
  Grad3 permGrad3[PSIZE];
  Grad4 permGrad4[PSIZE];

  /**
   * 2D SuperSimplex noise base.
   * Lookup table implementation inspired by DigitalShadow.
   */
  double noise2_Base(double xs, double ys);

  /**
   * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
   * Lookup table implementation inspired by DigitalShadow.
   * It was actually faster to narrow down the points in the loop itself,
   * than to build up the index with enough info to isolate 8 points.
   */
  double noise3_BCC(double xr, double yr, double zr);

  /**
   * 4D SuperSimplex noise base.
   * Using ultra-simple 4x4x4x4 lookup partitioning.
   * This isn't as elegant or SIMD/GPU/etc. portable as other approaches,
   * but it does compete performance-wise with optimized OpenSimplex1.
   */
  double noise4_Base(double xs, double ys, double zs, double ws);

  static int fastFloor(double x);

  static void initLatticePoints();
  static void initGradients();

  struct Initializer
  {
    Initializer();
  };
  static Initializer initializer;

public:
  explicit OpenSimplex2S(long seed = 0);

  /**
   * 2D SuperSimplex noise, standard lattice orientation.
   */
  double noise2(double x, double y);

  /**
   * 2D SuperSimplex noise, with Y pointing down the main diagonal.
   * Might be better for a 2D sandbox style game, where Y is vertical.
   * Probably slightly less optimal for heightmaps or continent maps.
   */
  double noise2_XBeforeY(double x, double y);

  /**
   * 3D Re-oriented 8-point BCC noise, classic orientation
   * Proper substitute for what 3D SuperSimplex would be,
   * in light of Forbidden Formulae.
   * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
   */
  double noise3_Classic(double x, double y, double z);

  /**
   * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Y).
   * Recommended for 3D terrain and time-varied animations.
   * The Z coordinate should always be the "different" coordinate in your use case.
   * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
   * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
   * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
   */
  double noise3_XYBeforeZ(double x, double y, double z);

  /**
   * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Z).
   * Recommended for 3D terrain and time-varied animations.
   * The Y coordinate should always be the "different" coordinate in your use case.
   * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
   * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
   * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
   */
  double noise3_XZBeforeY(double x, double y, double z);

  /**
   * 4D SuperSimplex noise, classic lattice orientation.
   */
  double noise4_Classic(double x, double y, double z, double w);

  /**
   * 4D SuperSimplex noise, with XY and ZW forming orthogonal triangular-based planes.
   * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
   * Recommended for noise(x, y, sin(time), cos(time)) trick.
   */
  double noise4_XYBeforeZW(double x, double y, double z, double w);

  /**
   * 4D SuperSimplex noise, with XZ and YW forming orthogonal triangular-based planes.
   * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
   */
  double noise4_XZBeforeYW(double x, double y, double z, double w);

  /**
   * 4D SuperSimplex noise, with XYZ oriented like noise3_Classic,
   * and W for an extra degree of freedom.
   * Recommended for time-varied animations which texture a 3D object (W=time)
   */
  double noise4_XYZBeforeW(double x, double y, double z, double w);
};
