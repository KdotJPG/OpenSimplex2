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

using System.Runtime.CompilerServices;

namespace Noise
{
    public class OpenSimplex2S
    {
        private const int PSIZE = 2048;
        private const int PMASK = 2047;

        private short[] perm;
        private Grad2[] permGrad2;
        private Grad3[] permGrad3;
        private Grad4[] permGrad4;

        public OpenSimplex2S(long seed)
        {
            perm = new short[PSIZE];
            permGrad2 = new Grad2[PSIZE];
            permGrad3 = new Grad3[PSIZE];
            permGrad4 = new Grad4[PSIZE];
            short[] source = new short[PSIZE];
            for (short i = 0; i < PSIZE; i++)
                source[i] = i;
            for (int i = PSIZE - 1; i >= 0; i--)
            {
                seed = seed * 6364136223846793005L + 1442695040888963407L;
                int r = (int)((seed + 31) % (i + 1));
                if (r < 0)
                    r += (i + 1);
                perm[i] = source[r];
                permGrad2[i] = GRADIENTS_2D[perm[i]];
                permGrad3[i] = GRADIENTS_3D[perm[i]];
                permGrad4[i] = GRADIENTS_4D[perm[i]];
                source[r] = source[i];
            }
        }

        /*
         * Noise Evaluators
         */

        /**
         * 2D SuperSimplex noise, standard lattice orientation.
         */
        public double Noise2(double x, double y)
        {

            // Get points for A2* lattice
            double s = 0.366025403784439 * (x + y);
            double xs = x + s, ys = y + s;

            return noise2_Base(xs, ys);
        }

        /**
         * 2D SuperSimplex noise, with Y pointing down the main diagonal.
         * Might be better for a 2D sandbox style game, where Y is vertical.
         * Probably slightly less optimal for heightmaps or continent maps.
         */
        public double Noise2_XBeforeY(double x, double y)
        {

            // Skew transform and rotation baked into one.
            double xx = x * 0.7071067811865476;
            double yy = y * 1.224744871380249;

            return noise2_Base(yy + xx, yy - xx);
        }

        /**
         * 2D SuperSimplex noise base.
         * Lookup table implementation inspired by DigitalShadow.
         */
        private double noise2_Base(double xs, double ys)
        {
            double value = 0;

            // Get base points and offsets
            int xsb = fastFloor(xs), ysb = fastFloor(ys);
            double xsi = xs - xsb, ysi = ys - ysb;

            // Index to point list
            int a = (int)(xsi + ysi);
            int index =
                (a << 2) |
                (int)(xsi - ysi / 2 + 1 - a / 2.0) << 3 |
                (int)(ysi - xsi / 2 + 1 - a / 2.0) << 4;

            double ssi = (xsi + ysi) * -0.211324865405187;
            double xi = xsi + ssi, yi = ysi + ssi;

            // Point contributions
            for (int i = 0; i < 4; i++)
            {
                LatticePoint2D c = LOOKUP_2D[index + i];

                double dx = xi + c.dx, dy = yi + c.dy;
                double attn = 2.0 / 3.0 - dx * dx - dy * dy;
                if (attn <= 0) continue;

                int pxm = (xsb + c.xsv) & PMASK, pym = (ysb + c.ysv) & PMASK;
                Grad2 grad = permGrad2[perm[pxm] ^ pym];
                double extrapolation = grad.dx * dx + grad.dy * dy;

                attn *= attn;
                value += attn * attn * extrapolation;
            }

            return value;
        }

        /**
         * 3D Re-oriented 8-point BCC noise, classic orientation
         * Proper substitute for what 3D SuperSimplex "should" be,
         * in light of Forbidden Formulae.
         * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
         */
        public double Noise3_Classic(double x, double y, double z)
        {

            // Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
            // If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
            // Orthonormal rotation. Not a skew transform.
            double r = (2.0 / 3.0) * (x + y + z);
            double xr = r - x, yr = r - y, zr = r - z;

            // Evaluate both lattices to form a BCC lattice.
            return noise3_BCC(xr, yr, zr);
        }

        /**
         * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Y).
         * Recommended for 3D terrain and time-varied animations.
         * The Z coordinate should always be the "different" coordinate in your use case.
         * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
         * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
         * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
         */
        public double Noise3_XYBeforeZ(double x, double y, double z)
        {

            // Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
            // Orthonormal rotation. Not a skew transform.
            double xy = x + y;
            double s2 = xy * -0.211324865405187;
            double zz = z * 0.577350269189626;
            double xr = x + s2 - zz, yr = y + s2 - zz;
            double zr = xy * 0.577350269189626 + zz;

            // Evaluate both lattices to form a BCC lattice.
            return noise3_BCC(xr, yr, zr);
        }

        /**
         * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Z).
         * Recommended for 3D terrain and time-varied animations.
         * The Y coordinate should always be the "different" coordinate in your use case.
         * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
         * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
         * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
         */
        public double Noise3_XZBeforeY(double x, double y, double z)
        {

            // Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
            // Orthonormal rotation. Not a skew transform.
            double xz = x + z;
            double s2 = xz * -0.211324865405187;
            double yy = y * 0.577350269189626;
            double xr = x + s2 - yy; double zr = z + s2 - yy;
            double yr = xz * 0.577350269189626 + yy;

            // Evaluate both lattices to form a BCC lattice.
            return noise3_BCC(xr, yr, zr);
        }

        /**
         * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
         * Lookup table implementation inspired by DigitalShadow.
         * It was actually faster to narrow down the points in the loop itself,
         * than to build up the index with enough info to isolate 8 points.
         */
        private double noise3_BCC(double xr, double yr, double zr)
        {

            // Get base and offsets inside cube of first lattice.
            int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr);
            double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

            // Identify which octant of the cube we're in. This determines which cell
            // in the other cubic lattice we're in, and also narrows down one point on each.
            int xht = (int)(xri + 0.5), yht = (int)(yri + 0.5), zht = (int)(zri + 0.5);
            int index = (xht << 0) | (yht << 1) | (zht << 2);

            // Point contributions
            double value = 0;
            LatticePoint3D c = LOOKUP_3D[index];
            while (c != null)
            {
                double dxr = xri + c.dxr, dyr = yri + c.dyr, dzr = zri + c.dzr;
                double attn = 0.75 - dxr * dxr - dyr * dyr - dzr * dzr;
                if (attn < 0)
                {
                    c = c.NextOnFailure;
                }
                else
                {
                    int pxm = (xrb + c.xrv) & PMASK, pym = (yrb + c.yrv) & PMASK, pzm = (zrb + c.zrv) & PMASK;
                    Grad3 grad = permGrad3[perm[perm[pxm] ^ pym] ^ pzm];
                    double extrapolation = grad.dx * dxr + grad.dy * dyr + grad.dz * dzr;

                    attn *= attn;
                    value += attn * attn * extrapolation;
                    c = c.NextOnSuccess;
                }
            }
            return value;
        }

        /**
		 * 4D SuperSimplex noise, classic lattice orientation.
		 */
        public double Noise4_Classic(double x, double y, double z, double w)
        {

            // Get points for A4 lattice
            double s = 0.309016994374947 * (x + y + z + w);
            double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

            return noise4_Base(xs, ys, zs, ws);
        }

        /**
		 * 4D SuperSimplex noise, with XY and ZW forming orthogonal triangular-based planes.
		 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
		 * Recommended for noise(x, y, sin(time), cos(time)) trick.
		 */
        public double Noise4_XYBeforeZW(double x, double y, double z, double w)
        {

            double s2 = (x + y) * -0.28522513987434876941 + (z + w) * 0.83897065470611435718;
            double t2 = (z + w) * 0.21939749883706435719 + (x + y) * -0.48214856493302476942;
            double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;

            return noise4_Base(xs, ys, zs, ws);
        }

        /**
		 * 4D SuperSimplex noise, with XZ and YW forming orthogonal triangular-based planes.
		 * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
		 */
        public double Noise4_XZBeforeYW(double x, double y, double z, double w)
        {

            double s2 = (x + z) * -0.28522513987434876941 + (y + w) * 0.83897065470611435718;
            double t2 = (y + w) * 0.21939749883706435719 + (x + z) * -0.48214856493302476942;
            double xs = x + s2, ys = y + t2, zs = z + s2, ws = w + t2;

            return noise4_Base(xs, ys, zs, ws);
        }

        /**
		 * 4D SuperSimplex noise, with XYZ oriented like noise3_Classic,
		 * and W for an extra degree of freedom.
		 * Recommended for time-varied animations which texture a 3D object (W=time)
		 */
        public double Noise4_XYZBeforeW(double x, double y, double z, double w)
        {

            double xyz = x + y + z;
            double ww = w * 1.118033988749894;
            double s2 = xyz * -0.16666666666666666 + ww;
            double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

            return noise4_Base(xs, ys, zs, ws);
        }

        /**
		 * 4D SuperSimplex noise base.
		 * Using ultra-simple 4x4x4x4 lookup partitioning.
		 * This isn't as elegant or SIMD/GPU/etc. portable as other approaches,
		 * but it does compete performance-wise with optimized OpenSimplex1.
		 */
        private double noise4_Base(double xs, double ys, double zs, double ws)
        {
            double value = 0;

            // Get base points and offsets
            int xsb = fastFloor(xs), ysb = fastFloor(ys), zsb = fastFloor(zs), wsb = fastFloor(ws);
            double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;

            // Unskewed offsets
            double ssi = (xsi + ysi + zsi + wsi) * -0.138196601125011;
            double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;

            int index = ((fastFloor(xs * 4) & 3) << 0)
                | ((fastFloor(ys * 4) & 3) << 2)
                | ((fastFloor(zs * 4) & 3) << 4)
                | ((fastFloor(ws * 4) & 3) << 6);

            // Point contributions
            foreach (LatticePoint4D c in LOOKUP_4D[index])
            {
                double dx = xi + c.dx, dy = yi + c.dy, dz = zi + c.dz, dw = wi + c.dw;
                double attn = 0.8 - dx * dx - dy * dy - dz * dz - dw * dw;
                if (attn > 0)
                {
                    attn *= attn;

                    int pxm = (xsb + c.xsv) & PMASK, pym = (ysb + c.ysv) & PMASK;
                    int pzm = (zsb + c.zsv) & PMASK, pwm = (wsb + c.wsv) & PMASK;
                    Grad4 grad = permGrad4[perm[perm[perm[pxm] ^ pym] ^ pzm] ^ pwm];
                    double extrapolation = grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

                    value += attn * attn * extrapolation;
                }
            }
            return value;
        }

        /*
         * Utility
         */

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int fastFloor(double x)
        {
            int xi = (int)x;
            return x < xi ? xi - 1 : xi;
        }

        /*
         * Lookup Tables & Gradients
         */

        private static LatticePoint2D[] LOOKUP_2D;
        private static LatticePoint3D[] LOOKUP_3D;
        private static LatticePoint4D[][] LOOKUP_4D;

        private const double N2 = 0.05481866495625118;
        private const double N3 = 0.2781926117527186;
        private const double N4 = 0.11127401889945551;
        private static Grad2[] GRADIENTS_2D;
        private static Grad3[] GRADIENTS_3D;
        private static Grad4[] GRADIENTS_4D;

        static OpenSimplex2S()
        {
            LOOKUP_2D = new LatticePoint2D[8 * 4];
            LOOKUP_3D = new LatticePoint3D[8];
            LOOKUP_4D = new LatticePoint4D[256][];

            for (int i = 0; i < 8; i++)
            {
                int i1, j1, i2, j2;
                if ((i & 1) == 0)
                {
                    if ((i & 2) == 0) { i1 = -1; j1 = 0; } else { i1 = 1; j1 = 0; }
                    if ((i & 4) == 0) { i2 = 0; j2 = -1; } else { i2 = 0; j2 = 1; }
                }
                else
                {
                    if ((i & 2) != 0) { i1 = 2; j1 = 1; } else { i1 = 0; j1 = 1; }
                    if ((i & 4) != 0) { i2 = 1; j2 = 2; } else { i2 = 1; j2 = 0; }
                }
                LOOKUP_2D[i * 4 + 0] = new LatticePoint2D(0, 0);
                LOOKUP_2D[i * 4 + 1] = new LatticePoint2D(1, 1);
                LOOKUP_2D[i * 4 + 2] = new LatticePoint2D(i1, j1);
                LOOKUP_2D[i * 4 + 3] = new LatticePoint2D(i2, j2);
            }

            for (int i = 0; i < 8; i++)
            {
                int i1, j1, k1, i2, j2, k2;
                i1 = (i >> 0) & 1; j1 = (i >> 1) & 1; k1 = (i >> 2) & 1;
                i2 = i1 ^ 1; j2 = j1 ^ 1; k2 = k1 ^ 1;

                // The two points within this octant, one from each of the two cubic half-lattices.
                LatticePoint3D c0 = new LatticePoint3D(i1, j1, k1, 0);
                LatticePoint3D c1 = new LatticePoint3D(i1 + i2, j1 + j2, k1 + k2, 1);

                // (1, 0, 0) vs (0, 1, 1) away from octant.
                LatticePoint3D c2 = new LatticePoint3D(i1 ^ 1, j1, k1, 0);
                LatticePoint3D c3 = new LatticePoint3D(i1, j1 ^ 1, k1 ^ 1, 0);

                // (1, 0, 0) vs (0, 1, 1) away from octant, on second half-lattice.
                LatticePoint3D c4 = new LatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + k2, 1);
                LatticePoint3D c5 = new LatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + (k2 ^ 1), 1);

                // (0, 1, 0) vs (1, 0, 1) away from octant.
                LatticePoint3D c6 = new LatticePoint3D(i1, j1 ^ 1, k1, 0);
                LatticePoint3D c7 = new LatticePoint3D(i1 ^ 1, j1, k1 ^ 1, 0);

                // (0, 1, 0) vs (1, 0, 1) away from octant, on second half-lattice.
                LatticePoint3D c8 = new LatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + k2, 1);
                LatticePoint3D c9 = new LatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + (k2 ^ 1), 1);

                // (0, 0, 1) vs (1, 1, 0) away from octant.
                LatticePoint3D cA = new LatticePoint3D(i1, j1, k1 ^ 1, 0);
                LatticePoint3D cB = new LatticePoint3D(i1 ^ 1, j1 ^ 1, k1, 0);

                // (0, 0, 1) vs (1, 1, 0) away from octant, on second half-lattice.
                LatticePoint3D cC = new LatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 ^ 1), 1);
                LatticePoint3D cD = new LatticePoint3D(i1 + (i2 ^ 1), j1 + (j2 ^ 1), k1 + k2, 1);

                // First two points are guaranteed.
                c0.NextOnFailure = c0.NextOnSuccess = c1;
                c1.NextOnFailure = c1.NextOnSuccess = c2;

                // If c2 is in range, then we know c3 and c4 are not.
                c2.NextOnFailure = c3; c2.NextOnSuccess = c5;
                c3.NextOnFailure = c4; c3.NextOnSuccess = c4;

                // If c4 is in range, then we know c5 is not.
                c4.NextOnFailure = c5; c4.NextOnSuccess = c6;
                c5.NextOnFailure = c5.NextOnSuccess = c6;

                // If c6 is in range, then we know c7 and c8 are not.
                c6.NextOnFailure = c7; c6.NextOnSuccess = c9;
                c7.NextOnFailure = c8; c7.NextOnSuccess = c8;

                // If c8 is in range, then we know c9 is not.
                c8.NextOnFailure = c9; c8.NextOnSuccess = cA;
                c9.NextOnFailure = c9.NextOnSuccess = cA;

                // If cA is in range, then we know cB and cC are not.
                cA.NextOnFailure = cB; cA.NextOnSuccess = cD;
                cB.NextOnFailure = cC; cB.NextOnSuccess = cC;

                // If cC is in range, then we know cD is not.
                cC.NextOnFailure = cD; cC.NextOnSuccess = null;
                cD.NextOnFailure = cD.NextOnSuccess = null;

                LOOKUP_3D[i] = c0;
            }

            int[][] lookup4DPregen = {
                new int[] { 0x15, 0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA },
                new int[] { 0x01, 0x05, 0x11, 0x15, 0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA },
                new int[] { 0x01, 0x15, 0x16, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA },
                new int[] { 0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA },
                new int[] { 0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA },
                new int[] { 0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA, 0xAB },
                new int[] { 0x04, 0x05, 0x14, 0x15, 0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA },
                new int[] { 0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA },
                new int[] { 0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA },
                new int[] { 0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x6A, 0x9A, 0xAA, 0xAB },
                new int[] { 0x04, 0x15, 0x19, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA, 0xAE },
                new int[] { 0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x5E, 0x6A, 0x9A, 0xAA, 0xAE },
                new int[] { 0x05, 0x15, 0x1A, 0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x5E, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x11, 0x15, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA },
                new int[] { 0x11, 0x15, 0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA },
                new int[] { 0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA, 0xAB },
                new int[] { 0x14, 0x15, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA6, 0xA9, 0xAA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x15, 0x16, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x6B, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x14, 0x15, 0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xAA },
                new int[] { 0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x6B, 0x9A, 0xAA, 0xAB },
                new int[] { 0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x19, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x9A, 0xAA, 0xAE },
                new int[] { 0x15, 0x1A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x10, 0x11, 0x14, 0x15, 0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA },
                new int[] { 0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA },
                new int[] { 0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x67, 0x6A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA },
                new int[] { 0x15, 0x16, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x6B, 0xA6, 0xAA, 0xAB },
                new int[] { 0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA },
                new int[] { 0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0xAA, 0xAB },
                new int[] { 0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x6D, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x19, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x6E, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0xAA, 0xAE },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x10, 0x15, 0x25, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA, 0xBA },
                new int[] { 0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0x76, 0xA6, 0xAA, 0xBA },
                new int[] { 0x11, 0x15, 0x26, 0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x67, 0x6A, 0x76, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB },
                new int[] { 0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x25, 0x55, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x25, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA6, 0xAA, 0xBA },
                new int[] { 0x15, 0x26, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB },
                new int[] { 0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0x79, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x25, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x25, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xAA, 0xBA },
                new int[] { 0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB },
                new int[] { 0x14, 0x15, 0x29, 0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0x6D, 0x79, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE },
                new int[] { 0x15, 0x29, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE },
                new int[] { 0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x7A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF },
                new int[] { 0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x41, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA },
                new int[] { 0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA },
                new int[] { 0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x44, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x45, 0x46, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB },
                new int[] { 0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA },
                new int[] { 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x9A, 0x9B, 0xAA, 0xAB },
                new int[] { 0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x45, 0x49, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE },
                new int[] { 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x99, 0x9A, 0x9E, 0xAA, 0xAE },
                new int[] { 0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x51, 0x55, 0x56, 0x59, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB },
                new int[] { 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x56, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE },
                new int[] { 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA },
                new int[] { 0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xA7, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA },
                new int[] { 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA },
                new int[] { 0x15, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAD, 0xAE },
                new int[] { 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x51, 0x55, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA },
                new int[] { 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA },
                new int[] { 0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x6A, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB },
                new int[] { 0x54, 0x55, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA },
                new int[] { 0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB },
                new int[] { 0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA },
                new int[] { 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE },
                new int[] { 0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE },
                new int[] { 0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF },
                new int[] { 0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA },
                new int[] { 0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA },
                new int[] { 0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA },
                new int[] { 0x45, 0x46, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB },
                new int[] { 0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA },
                new int[] { 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xAA, 0xAB },
                new int[] { 0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x95, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE },
                new int[] { 0x45, 0x49, 0x55, 0x59, 0x5A, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE },
                new int[] { 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xAA, 0xAE },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA },
                new int[] { 0x51, 0x52, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xEA },
                new int[] { 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA },
                new int[] { 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xEA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x58, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE },
                new int[] { 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA },
                new int[] { 0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB },
                new int[] { 0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE },
                new int[] { 0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xAF },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA },
                new int[] { 0x51, 0x55, 0x61, 0x65, 0x66, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA },
                new int[] { 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB },
                new int[] { 0x54, 0x55, 0x64, 0x65, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA },
                new int[] { 0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB },
                new int[] { 0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA },
                new int[] { 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xBB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE },
                new int[] { 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA },
                new int[] { 0x40, 0x45, 0x51, 0x54, 0x55, 0x85, 0x91, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA, 0xEA },
                new int[] { 0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xD6, 0xEA },
                new int[] { 0x41, 0x45, 0x51, 0x55, 0x56, 0x86, 0x92, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB, 0xD6, 0xEA, 0xEB },
                new int[] { 0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xEA },
                new int[] { 0x45, 0x55, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xDA, 0xEA },
                new int[] { 0x45, 0x55, 0x56, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xDA, 0xEA },
                new int[] { 0x45, 0x55, 0x56, 0x86, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB },
                new int[] { 0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xD9, 0xEA },
                new int[] { 0x45, 0x55, 0x59, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xDA, 0xEA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xDA, 0xEA },
                new int[] { 0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB },
                new int[] { 0x44, 0x45, 0x54, 0x55, 0x59, 0x89, 0x95, 0x98, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE, 0xD9, 0xEA, 0xEE },
                new int[] { 0x45, 0x55, 0x59, 0x89, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE },
                new int[] { 0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xDA, 0xEA, 0xEF },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x51, 0x55, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA },
                new int[] { 0x51, 0x55, 0x56, 0x91, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA },
                new int[] { 0x51, 0x55, 0x56, 0x92, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB },
                new int[] { 0x54, 0x55, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA },
                new int[] { 0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB },
                new int[] { 0x54, 0x55, 0x59, 0x94, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA },
                new int[] { 0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB },
                new int[] { 0x54, 0x55, 0x59, 0x95, 0x98, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE },
                new int[] { 0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE },
                new int[] { 0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xEA, 0xEF },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xE5, 0xEA },
                new int[] { 0x51, 0x55, 0x65, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA },
                new int[] { 0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB },
                new int[] { 0x54, 0x55, 0x65, 0x94, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA },
                new int[] { 0x55, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x94, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x55, 0x56, 0x59, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA },
                new int[] { 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xEA, 0xEB },
                new int[] { 0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE },
                new int[] { 0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE },
                new int[] { 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xEA },
                new int[] { 0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA1, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA, 0xE5, 0xEA, 0xFA },
                new int[] { 0x51, 0x55, 0x65, 0x95, 0xA1, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA },
                new int[] { 0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA },
                new int[] { 0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xE6, 0xEA, 0xFB },
                new int[] { 0x54, 0x55, 0x65, 0x95, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA },
                new int[] { 0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA },
                new int[] { 0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA },
                new int[] { 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xEA, 0xFB },
                new int[] { 0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA },
                new int[] { 0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA },
                new int[] { 0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA },
                new int[] { 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xEA },
                new int[] { 0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xE9, 0xEA, 0xFE },
                new int[] { 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA, 0xFE },
                new int[] { 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA },
                new int[] { 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA, 0xEA },
            };
            LatticePoint4D[] latticePoints = new LatticePoint4D[256];
            for (int i = 0; i < 256; i++)
            {
                int cx = ((i >> 0) & 3) - 1;
                int cy = ((i >> 2) & 3) - 1;
                int cz = ((i >> 4) & 3) - 1;
                int cw = ((i >> 6) & 3) - 1;
                latticePoints[i] = new LatticePoint4D(cx, cy, cz, cw);
            }
            for (int i = 0; i < 256; i++)
            {
                LOOKUP_4D[i] = new LatticePoint4D[lookup4DPregen[i].Length];
                for (int j = 0; j < lookup4DPregen[i].Length; j++)
                {
                    LOOKUP_4D[i][j] = latticePoints[lookup4DPregen[i][j]];
                }
            }

            GRADIENTS_2D = new Grad2[PSIZE];
            Grad2[] grad2 = {
                new Grad2( 0.130526192220052,  0.99144486137381),
                new Grad2( 0.38268343236509,   0.923879532511287),
                new Grad2( 0.608761429008721,  0.793353340291235),
                new Grad2( 0.793353340291235,  0.608761429008721),
                new Grad2( 0.923879532511287,  0.38268343236509),
                new Grad2( 0.99144486137381,   0.130526192220051),
                new Grad2( 0.99144486137381,  -0.130526192220051),
                new Grad2( 0.923879532511287, -0.38268343236509),
                new Grad2( 0.793353340291235, -0.60876142900872),
                new Grad2( 0.608761429008721, -0.793353340291235),
                new Grad2( 0.38268343236509,  -0.923879532511287),
                new Grad2( 0.130526192220052, -0.99144486137381),
                new Grad2(-0.130526192220052, -0.99144486137381),
                new Grad2(-0.38268343236509,  -0.923879532511287),
                new Grad2(-0.608761429008721, -0.793353340291235),
                new Grad2(-0.793353340291235, -0.608761429008721),
                new Grad2(-0.923879532511287, -0.38268343236509),
                new Grad2(-0.99144486137381,  -0.130526192220052),
                new Grad2(-0.99144486137381,   0.130526192220051),
                new Grad2(-0.923879532511287,  0.38268343236509),
                new Grad2(-0.793353340291235,  0.608761429008721),
                new Grad2(-0.608761429008721,  0.793353340291235),
                new Grad2(-0.38268343236509,   0.923879532511287),
                new Grad2(-0.130526192220052,  0.99144486137381)
            };
            for (int i = 0; i < grad2.Length; i++)
            {
                grad2[i].dx /= N2; grad2[i].dy /= N2;
            }
            for (int i = 0; i < PSIZE; i++)
            {
                GRADIENTS_2D[i] = grad2[i % grad2.Length];
            }

            GRADIENTS_3D = new Grad3[PSIZE];
            Grad3[] grad3 = {
                new Grad3(-2.22474487139,      -2.22474487139,      -1.0),
                new Grad3(-2.22474487139,      -2.22474487139,       1.0),
                new Grad3(-3.0862664687972017, -1.1721513422464978,  0.0),
                new Grad3(-1.1721513422464978, -3.0862664687972017,  0.0),
                new Grad3(-2.22474487139,      -1.0,                -2.22474487139),
                new Grad3(-2.22474487139,       1.0,                -2.22474487139),
                new Grad3(-1.1721513422464978,  0.0,                -3.0862664687972017),
                new Grad3(-3.0862664687972017,  0.0,                -1.1721513422464978),
                new Grad3(-2.22474487139,      -1.0,                 2.22474487139),
                new Grad3(-2.22474487139,       1.0,                 2.22474487139),
                new Grad3(-3.0862664687972017,  0.0,                 1.1721513422464978),
                new Grad3(-1.1721513422464978,  0.0,                 3.0862664687972017),
                new Grad3(-2.22474487139,       2.22474487139,      -1.0),
                new Grad3(-2.22474487139,       2.22474487139,       1.0),
                new Grad3(-1.1721513422464978,  3.0862664687972017,  0.0),
                new Grad3(-3.0862664687972017,  1.1721513422464978,  0.0),
                new Grad3(-1.0,                -2.22474487139,      -2.22474487139),
                new Grad3( 1.0,                -2.22474487139,      -2.22474487139),
                new Grad3( 0.0,                -3.0862664687972017, -1.1721513422464978),
                new Grad3( 0.0,                -1.1721513422464978, -3.0862664687972017),
                new Grad3(-1.0,                -2.22474487139,       2.22474487139),
                new Grad3( 1.0,                -2.22474487139,       2.22474487139),
                new Grad3( 0.0,                -1.1721513422464978,  3.0862664687972017),
                new Grad3( 0.0,                -3.0862664687972017,  1.1721513422464978),
                new Grad3(-1.0,                 2.22474487139,      -2.22474487139),
                new Grad3( 1.0,                 2.22474487139,      -2.22474487139),
                new Grad3( 0.0,                 1.1721513422464978, -3.0862664687972017),
                new Grad3( 0.0,                 3.0862664687972017, -1.1721513422464978),
                new Grad3(-1.0,                 2.22474487139,       2.22474487139),
                new Grad3( 1.0,                 2.22474487139,       2.22474487139),
                new Grad3( 0.0,                 3.0862664687972017,  1.1721513422464978),
                new Grad3( 0.0,                 1.1721513422464978,  3.0862664687972017),
                new Grad3( 2.22474487139,      -2.22474487139,      -1.0),
                new Grad3( 2.22474487139,      -2.22474487139,       1.0),
                new Grad3( 1.1721513422464978, -3.0862664687972017,  0.0),
                new Grad3( 3.0862664687972017, -1.1721513422464978,  0.0),
                new Grad3( 2.22474487139,      -1.0,                -2.22474487139),
                new Grad3( 2.22474487139,       1.0,                -2.22474487139),
                new Grad3( 3.0862664687972017,  0.0,                -1.1721513422464978),
                new Grad3( 1.1721513422464978,  0.0,                -3.0862664687972017),
                new Grad3( 2.22474487139,      -1.0,                 2.22474487139),
                new Grad3( 2.22474487139,       1.0,                 2.22474487139),
                new Grad3( 1.1721513422464978,  0.0,                 3.0862664687972017),
                new Grad3( 3.0862664687972017,  0.0,                 1.1721513422464978),
                new Grad3( 2.22474487139,       2.22474487139,      -1.0),
                new Grad3( 2.22474487139,       2.22474487139,       1.0),
                new Grad3( 3.0862664687972017,  1.1721513422464978,  0.0),
                new Grad3( 1.1721513422464978,  3.0862664687972017,  0.0)
            };
            for (int i = 0; i < grad3.Length; i++)
            {
                grad3[i].dx /= N3; grad3[i].dy /= N3; grad3[i].dz /= N3;
            }
            for (int i = 0; i < PSIZE; i++)
            {
                GRADIENTS_3D[i] = grad3[i % grad3.Length];
            }

            GRADIENTS_4D = new Grad4[PSIZE];
            Grad4[] grad4 = {
                new Grad4(-0.753341017856078,    -0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624),
                new Grad4(-0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301,    0.12128480194602098),
                new Grad4(-0.7821684431180708,   -0.4321472685365301,    0.12128480194602098,  -0.4321472685365301),
                new Grad4(-0.7821684431180708,    0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301),
                new Grad4(-0.8586508742123365,   -0.508629699630796,     0.044802370851755174,  0.044802370851755174),
                new Grad4(-0.8586508742123365,    0.044802370851755174, -0.508629699630796,     0.044802370851755174),
                new Grad4(-0.8586508742123365,    0.044802370851755174,  0.044802370851755174, -0.508629699630796),
                new Grad4(-0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842),
                new Grad4(-0.37968289875261624,  -0.753341017856078,    -0.37968289875261624,  -0.37968289875261624),
                new Grad4(-0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301,    0.12128480194602098),
                new Grad4(-0.4321472685365301,   -0.7821684431180708,    0.12128480194602098,  -0.4321472685365301),
                new Grad4( 0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301),
                new Grad4(-0.508629699630796,    -0.8586508742123365,    0.044802370851755174,  0.044802370851755174),
                new Grad4( 0.044802370851755174, -0.8586508742123365,   -0.508629699630796,     0.044802370851755174),
                new Grad4( 0.044802370851755174, -0.8586508742123365,    0.044802370851755174, -0.508629699630796),
                new Grad4(-0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842),
                new Grad4(-0.37968289875261624,  -0.37968289875261624,  -0.753341017856078,    -0.37968289875261624),
                new Grad4(-0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708,    0.12128480194602098),
                new Grad4(-0.4321472685365301,    0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301),
                new Grad4( 0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301),
                new Grad4(-0.508629699630796,     0.044802370851755174, -0.8586508742123365,    0.044802370851755174),
                new Grad4( 0.044802370851755174, -0.508629699630796,    -0.8586508742123365,    0.044802370851755174),
                new Grad4( 0.044802370851755174,  0.044802370851755174, -0.8586508742123365,   -0.508629699630796),
                new Grad4(-0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842),
                new Grad4(-0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624,  -0.753341017856078),
                new Grad4(-0.4321472685365301,   -0.4321472685365301,    0.12128480194602098,  -0.7821684431180708),
                new Grad4(-0.4321472685365301,    0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708),
                new Grad4( 0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708),
                new Grad4(-0.508629699630796,     0.044802370851755174,  0.044802370851755174, -0.8586508742123365),
                new Grad4( 0.044802370851755174, -0.508629699630796,     0.044802370851755174, -0.8586508742123365),
                new Grad4( 0.044802370851755174,  0.044802370851755174, -0.508629699630796,    -0.8586508742123365),
                new Grad4(-0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062),
                new Grad4(-0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537,    0.5794684678643381),
                new Grad4(-0.7504883828755602,   -0.4004672082940195,    0.15296486218853164,   0.5029860367700724),
                new Grad4(-0.7504883828755602,    0.15296486218853164,  -0.4004672082940195,    0.5029860367700724),
                new Grad4(-0.8828161875373585,    0.08164729285680945,   0.08164729285680945,   0.4553054119602712),
                new Grad4(-0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945,   0.8828161875373585),
                new Grad4(-0.5029860367700724,   -0.15296486218853164,   0.4004672082940195,    0.7504883828755602),
                new Grad4(-0.5029860367700724,    0.4004672082940195,   -0.15296486218853164,   0.7504883828755602),
                new Grad4(-0.5794684678643381,    0.3239847771997537,    0.3239847771997537,    0.6740059517812944),
                new Grad4(-0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537,    0.5794684678643381),
                new Grad4(-0.4004672082940195,   -0.7504883828755602,    0.15296486218853164,   0.5029860367700724),
                new Grad4( 0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195,    0.5029860367700724),
                new Grad4( 0.08164729285680945,  -0.8828161875373585,    0.08164729285680945,   0.4553054119602712),
                new Grad4(-0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945,   0.8828161875373585),
                new Grad4(-0.15296486218853164,  -0.5029860367700724,    0.4004672082940195,    0.7504883828755602),
                new Grad4( 0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164,   0.7504883828755602),
                new Grad4( 0.3239847771997537,   -0.5794684678643381,    0.3239847771997537,    0.6740059517812944),
                new Grad4(-0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944,    0.5794684678643381),
                new Grad4(-0.4004672082940195,    0.15296486218853164,  -0.7504883828755602,    0.5029860367700724),
                new Grad4( 0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602,    0.5029860367700724),
                new Grad4( 0.08164729285680945,   0.08164729285680945,  -0.8828161875373585,    0.4553054119602712),
                new Grad4(-0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712,    0.8828161875373585),
                new Grad4(-0.15296486218853164,   0.4004672082940195,   -0.5029860367700724,    0.7504883828755602),
                new Grad4( 0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724,    0.7504883828755602),
                new Grad4( 0.3239847771997537,    0.3239847771997537,   -0.5794684678643381,    0.6740059517812944),
                new Grad4(-0.6740059517812944,   -0.3239847771997537,    0.5794684678643381,   -0.3239847771997537),
                new Grad4(-0.7504883828755602,   -0.4004672082940195,    0.5029860367700724,    0.15296486218853164),
                new Grad4(-0.7504883828755602,    0.15296486218853164,   0.5029860367700724,   -0.4004672082940195),
                new Grad4(-0.8828161875373585,    0.08164729285680945,   0.4553054119602712,    0.08164729285680945),
                new Grad4(-0.4553054119602712,   -0.08164729285680945,   0.8828161875373585,   -0.08164729285680945),
                new Grad4(-0.5029860367700724,   -0.15296486218853164,   0.7504883828755602,    0.4004672082940195),
                new Grad4(-0.5029860367700724,    0.4004672082940195,    0.7504883828755602,   -0.15296486218853164),
                new Grad4(-0.5794684678643381,    0.3239847771997537,    0.6740059517812944,    0.3239847771997537),
                new Grad4(-0.3239847771997537,   -0.6740059517812944,    0.5794684678643381,   -0.3239847771997537),
                new Grad4(-0.4004672082940195,   -0.7504883828755602,    0.5029860367700724,    0.15296486218853164),
                new Grad4( 0.15296486218853164,  -0.7504883828755602,    0.5029860367700724,   -0.4004672082940195),
                new Grad4( 0.08164729285680945,  -0.8828161875373585,    0.4553054119602712,    0.08164729285680945),
                new Grad4(-0.08164729285680945,  -0.4553054119602712,    0.8828161875373585,   -0.08164729285680945),
                new Grad4(-0.15296486218853164,  -0.5029860367700724,    0.7504883828755602,    0.4004672082940195),
                new Grad4( 0.4004672082940195,   -0.5029860367700724,    0.7504883828755602,   -0.15296486218853164),
                new Grad4( 0.3239847771997537,   -0.5794684678643381,    0.6740059517812944,    0.3239847771997537),
                new Grad4(-0.3239847771997537,   -0.3239847771997537,    0.5794684678643381,   -0.6740059517812944),
                new Grad4(-0.4004672082940195,    0.15296486218853164,   0.5029860367700724,   -0.7504883828755602),
                new Grad4( 0.15296486218853164,  -0.4004672082940195,    0.5029860367700724,   -0.7504883828755602),
                new Grad4( 0.08164729285680945,   0.08164729285680945,   0.4553054119602712,   -0.8828161875373585),
                new Grad4(-0.08164729285680945,  -0.08164729285680945,   0.8828161875373585,   -0.4553054119602712),
                new Grad4(-0.15296486218853164,   0.4004672082940195,    0.7504883828755602,   -0.5029860367700724),
                new Grad4( 0.4004672082940195,   -0.15296486218853164,   0.7504883828755602,   -0.5029860367700724),
                new Grad4( 0.3239847771997537,    0.3239847771997537,    0.6740059517812944,   -0.5794684678643381),
                new Grad4(-0.6740059517812944,    0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537),
                new Grad4(-0.7504883828755602,    0.5029860367700724,   -0.4004672082940195,    0.15296486218853164),
                new Grad4(-0.7504883828755602,    0.5029860367700724,    0.15296486218853164,  -0.4004672082940195),
                new Grad4(-0.8828161875373585,    0.4553054119602712,    0.08164729285680945,   0.08164729285680945),
                new Grad4(-0.4553054119602712,    0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945),
                new Grad4(-0.5029860367700724,    0.7504883828755602,   -0.15296486218853164,   0.4004672082940195),
                new Grad4(-0.5029860367700724,    0.7504883828755602,    0.4004672082940195,   -0.15296486218853164),
                new Grad4(-0.5794684678643381,    0.6740059517812944,    0.3239847771997537,    0.3239847771997537),
                new Grad4(-0.3239847771997537,    0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537),
                new Grad4(-0.4004672082940195,    0.5029860367700724,   -0.7504883828755602,    0.15296486218853164),
                new Grad4( 0.15296486218853164,   0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195),
                new Grad4( 0.08164729285680945,   0.4553054119602712,   -0.8828161875373585,    0.08164729285680945),
                new Grad4(-0.08164729285680945,   0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945),
                new Grad4(-0.15296486218853164,   0.7504883828755602,   -0.5029860367700724,    0.4004672082940195),
                new Grad4( 0.4004672082940195,    0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164),
                new Grad4( 0.3239847771997537,    0.6740059517812944,   -0.5794684678643381,    0.3239847771997537),
                new Grad4(-0.3239847771997537,    0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944),
                new Grad4(-0.4004672082940195,    0.5029860367700724,    0.15296486218853164,  -0.7504883828755602),
                new Grad4( 0.15296486218853164,   0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602),
                new Grad4( 0.08164729285680945,   0.4553054119602712,    0.08164729285680945,  -0.8828161875373585),
                new Grad4(-0.08164729285680945,   0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712),
                new Grad4(-0.15296486218853164,   0.7504883828755602,    0.4004672082940195,   -0.5029860367700724),
                new Grad4( 0.4004672082940195,    0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724),
                new Grad4( 0.3239847771997537,    0.6740059517812944,    0.3239847771997537,   -0.5794684678643381),
                new Grad4( 0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537),
                new Grad4( 0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195,    0.15296486218853164),
                new Grad4( 0.5029860367700724,   -0.7504883828755602,    0.15296486218853164,  -0.4004672082940195),
                new Grad4( 0.4553054119602712,   -0.8828161875373585,    0.08164729285680945,   0.08164729285680945),
                new Grad4( 0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945),
                new Grad4( 0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164,   0.4004672082940195),
                new Grad4( 0.7504883828755602,   -0.5029860367700724,    0.4004672082940195,   -0.15296486218853164),
                new Grad4( 0.6740059517812944,   -0.5794684678643381,    0.3239847771997537,    0.3239847771997537),
                new Grad4( 0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537),
                new Grad4( 0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602,    0.15296486218853164),
                new Grad4( 0.5029860367700724,    0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195),
                new Grad4( 0.4553054119602712,    0.08164729285680945,  -0.8828161875373585,    0.08164729285680945),
                new Grad4( 0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945),
                new Grad4( 0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724,    0.4004672082940195),
                new Grad4( 0.7504883828755602,    0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164),
                new Grad4( 0.6740059517812944,    0.3239847771997537,   -0.5794684678643381,    0.3239847771997537),
                new Grad4( 0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944),
                new Grad4( 0.5029860367700724,   -0.4004672082940195,    0.15296486218853164,  -0.7504883828755602),
                new Grad4( 0.5029860367700724,    0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602),
                new Grad4( 0.4553054119602712,    0.08164729285680945,   0.08164729285680945,  -0.8828161875373585),
                new Grad4( 0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712),
                new Grad4( 0.7504883828755602,   -0.15296486218853164,   0.4004672082940195,   -0.5029860367700724),
                new Grad4( 0.7504883828755602,    0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724),
                new Grad4( 0.6740059517812944,    0.3239847771997537,    0.3239847771997537,   -0.5794684678643381),
                new Grad4( 0.03381941603233842,   0.03381941603233842,   0.03381941603233842,   0.9982828964265062),
                new Grad4(-0.044802370851755174, -0.044802370851755174,  0.508629699630796,     0.8586508742123365),
                new Grad4(-0.044802370851755174,  0.508629699630796,    -0.044802370851755174,  0.8586508742123365),
                new Grad4(-0.12128480194602098,   0.4321472685365301,    0.4321472685365301,    0.7821684431180708),
                new Grad4( 0.508629699630796,    -0.044802370851755174, -0.044802370851755174,  0.8586508742123365),
                new Grad4( 0.4321472685365301,   -0.12128480194602098,   0.4321472685365301,    0.7821684431180708),
                new Grad4( 0.4321472685365301,    0.4321472685365301,   -0.12128480194602098,   0.7821684431180708),
                new Grad4( 0.37968289875261624,   0.37968289875261624,   0.37968289875261624,   0.753341017856078),
                new Grad4( 0.03381941603233842,   0.03381941603233842,   0.9982828964265062,    0.03381941603233842),
                new Grad4(-0.044802370851755174,  0.044802370851755174,  0.8586508742123365,    0.508629699630796),
                new Grad4(-0.044802370851755174,  0.508629699630796,     0.8586508742123365,   -0.044802370851755174),
                new Grad4(-0.12128480194602098,   0.4321472685365301,    0.7821684431180708,    0.4321472685365301),
                new Grad4( 0.508629699630796,    -0.044802370851755174,  0.8586508742123365,   -0.044802370851755174),
                new Grad4( 0.4321472685365301,   -0.12128480194602098,   0.7821684431180708,    0.4321472685365301),
                new Grad4( 0.4321472685365301,    0.4321472685365301,    0.7821684431180708,   -0.12128480194602098),
                new Grad4( 0.37968289875261624,   0.37968289875261624,   0.753341017856078,     0.37968289875261624),
                new Grad4( 0.03381941603233842,   0.9982828964265062,    0.03381941603233842,   0.03381941603233842),
                new Grad4(-0.044802370851755174,  0.8586508742123365,   -0.044802370851755174,  0.508629699630796),
                new Grad4(-0.044802370851755174,  0.8586508742123365,    0.508629699630796,    -0.044802370851755174),
                new Grad4(-0.12128480194602098,   0.7821684431180708,    0.4321472685365301,    0.4321472685365301),
                new Grad4( 0.508629699630796,     0.8586508742123365,   -0.044802370851755174, -0.044802370851755174),
                new Grad4( 0.4321472685365301,    0.7821684431180708,   -0.12128480194602098,   0.4321472685365301),
                new Grad4( 0.4321472685365301,    0.7821684431180708,    0.4321472685365301,   -0.12128480194602098),
                new Grad4( 0.37968289875261624,   0.753341017856078,     0.37968289875261624,   0.37968289875261624),
                new Grad4( 0.9982828964265062,    0.03381941603233842,   0.03381941603233842,   0.03381941603233842),
                new Grad4( 0.8586508742123365,   -0.044802370851755174, -0.044802370851755174,  0.508629699630796),
                new Grad4( 0.8586508742123365,   -0.044802370851755174,  0.508629699630796,    -0.044802370851755174),
                new Grad4( 0.7821684431180708,   -0.12128480194602098,   0.4321472685365301,    0.4321472685365301),
                new Grad4( 0.8586508742123365,    0.508629699630796,    -0.044802370851755174, -0.044802370851755174),
                new Grad4( 0.7821684431180708,    0.4321472685365301,   -0.12128480194602098,   0.4321472685365301),
                new Grad4( 0.7821684431180708,    0.4321472685365301,    0.4321472685365301,   -0.12128480194602098),
                new Grad4( 0.753341017856078,     0.37968289875261624,   0.37968289875261624,   0.37968289875261624)
            };
            for (int i = 0; i < grad4.Length; i++)
            {
                grad4[i].dx /= N4; grad4[i].dy /= N4; grad4[i].dz /= N4; grad4[i].dw /= N4;
            }
            for (int i = 0; i < PSIZE; i++)
            {
                GRADIENTS_4D[i] = grad4[i % grad4.Length];
            }
        }

        private struct LatticePoint2D
        {
            public int xsv, ysv;
            public double dx, dy;
            public LatticePoint2D(int xsv, int ysv)
            {
                this.xsv = xsv; this.ysv = ysv;
                double ssv = (xsv + ysv) * -0.211324865405187;
                this.dx = -xsv - ssv;
                this.dy = -ysv - ssv;
            }
        }

        private class LatticePoint3D
        {
            public double dxr, dyr, dzr;
            public int xrv, yrv, zrv;
            public LatticePoint3D NextOnFailure, NextOnSuccess;
            public LatticePoint3D(int xrv, int yrv, int zrv, int lattice)
            {
                this.dxr = -xrv + lattice * 0.5; this.dyr = -yrv + lattice * 0.5; this.dzr = -zrv + lattice * 0.5;
                this.xrv = xrv + lattice * 1024; this.yrv = yrv + lattice * 1024; this.zrv = zrv + lattice * 1024;
            }
        }

        private class LatticePoint4D
        {
            public int xsv, ysv, zsv, wsv;
            public double dx, dy, dz, dw;
            public LatticePoint4D(int xsv, int ysv, int zsv, int wsv)
            {
                this.xsv = xsv; this.ysv = ysv; this.zsv = zsv; this.wsv = wsv;
                double ssv = (xsv + ysv + zsv + wsv) * -0.138196601125011;
                this.dx = -xsv - ssv;
                this.dy = -ysv - ssv;
                this.dz = -zsv - ssv;
                this.dw = -wsv - ssv;
            }
        }

        private struct Grad2
        {
            public double dx, dy;
            public Grad2(double dx, double dy)
            {
                this.dx = dx; this.dy = dy;
            }
        }

        private struct Grad3
        {
            public double dx, dy, dz;
            public Grad3(double dx, double dy, double dz)
            {
                this.dx = dx; this.dy = dy; this.dz = dz;
            }
        }

        private struct Grad4
        {
            public double dx, dy, dz, dw;
            public Grad4(double dx, double dy, double dz, double dw)
            {
                this.dx = dx; this.dy = dy; this.dz = dz; this.dw = dw;
            }
        }
    }
}
