/**
 * K.jpg's OpenSimplex 2, faster variant ("Fast Simplex-Style Noise")
 *
 * - 2D is standard simplex implemented using a lookup table.
 * - 3D is "Re-oriented 4-point BCC noise" which constructs an
 *   isomorphic BCC lattice in a much different way than usual.
 *
 * Multiple versions of each function are provided. See the
 * documentation above each, for more info.
 */

using System.Runtime.CompilerServices;

namespace Noise
{
    public class OpenSimplex2F
    {
        private const int PSIZE = 2048;
        private const int PMASK = 2047;

        private short[] perm;
        private Grad2[] permGrad2;
        private Grad3[] permGrad3;

        public OpenSimplex2F(long seed)
        {
            perm = new short[PSIZE];
            permGrad2 = new Grad2[PSIZE];
            permGrad3 = new Grad3[PSIZE];
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
                source[r] = source[i];
            }
        }

        /*
         * Noise Evaluators
         */

        /**
         * 2D Simplex noise, standard lattice orientation.
         */
        public double Noise2(double x, double y)
        {

            // Get points for A2* lattice
            double s = 0.366025403784439 * (x + y);
            double xs = x + s, ys = y + s;

            return noise2_Base(xs, ys);
        }

        /**
         * 2D Simplex noise, with Y pointing down the main diagonal.
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
         * 2D Simplex noise base.
         * Lookup table implementation inspired by DigitalShadow.
         */
        private double noise2_Base(double xs, double ys)
        {
            double value = 0;

            // Get base points and offsets
            int xsb = fastFloor(xs), ysb = fastFloor(ys);
            double xsi = xs - xsb, ysi = ys - ysb;

            // Index to point list
            int index = (int)((ysi - xsi) / 2 + 1) * 3;

            double ssi = (xsi + ysi) * -0.211324865405187;
            double xi = xsi + ssi, yi = ysi + ssi;

            // Point contributions
            for (int i = 0; i < 3; i++)
            {
                LatticePoint2D c = LOOKUP_2D[index + i];

                double dx = xi + c.dx, dy = yi + c.dy;
                double attn = 0.5 - dx * dx - dy * dy;
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
         * 3D Re-oriented 4-point BCC noise, classic orientation.
         * Proper substitute for 3D Simplex in light of Forbidden Formulae.
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
         * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Y).
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
         * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Z).
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
         * than to build up the index with enough info to isolate 4 points.
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
                double attn = 0.5 - dxr * dxr - dyr * dyr - dzr * dzr;
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

        private const double N2 = 0.01001634121365712;
        private const double N3 = 0.030485933181293584;
        private static Grad2[] GRADIENTS_2D;
        private static Grad3[] GRADIENTS_3D;

        static OpenSimplex2F() {
            LOOKUP_2D = new LatticePoint2D[2 * 3];
            LOOKUP_3D = new LatticePoint3D[8];
        
            for (int i = 0; i < 2; i++) {
                int i1, j1;
                if ((i & 1) == 0) { i1 = 1; j1 = 0; }
                else { i1 = 0; j1 = 1; }
                LOOKUP_2D[i * 3 + 0] = new LatticePoint2D(0, 0);
                LOOKUP_2D[i * 3 + 1] = new LatticePoint2D(1, 1);
                LOOKUP_2D[i * 3 + 2] = new LatticePoint2D(i1, j1);
            }
        
            for (int i = 0; i < 8; i++) {
                int i1, j1, k1, i2, j2, k2;
                i1 = (i >> 0) & 1; j1 = (i >> 1) & 1; k1 = (i >> 2) & 1;
                i2 = i1 ^ 1; j2 = j1 ^ 1; k2 = k1 ^ 1;
            
                // The two points within this octant, one from each of the two cubic half-lattices.
                LatticePoint3D c0 = new LatticePoint3D(i1, j1, k1, 0);
                LatticePoint3D c1 = new LatticePoint3D(i1 + i2, j1 + j2, k1 + k2, 1);

                // Each single step away on the first half-lattice.
                LatticePoint3D c2 = new LatticePoint3D(i1 ^ 1, j1, k1, 0);
                LatticePoint3D c3 = new LatticePoint3D(i1, j1 ^ 1, k1, 0);
                LatticePoint3D c4 = new LatticePoint3D(i1, j1, k1 ^ 1, 0);

                // Each single step away on the second half-lattice.
                LatticePoint3D c5 = new LatticePoint3D(i1 + (i2 ^ 1), j1 + j2, k1 + k2, 1);
                LatticePoint3D c6 = new LatticePoint3D(i1 + i2, j1 + (j2 ^ 1), k1 + k2, 1);
                LatticePoint3D c7 = new LatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 ^ 1), 1);

                // First two are guaranteed.
                c0.NextOnFailure = c0.NextOnSuccess = c1;
                c1.NextOnFailure = c1.NextOnSuccess = c2;
            
                // Once we find one on the first half-lattice, the rest are out.
                // In addition, knowing c2 rules out c5.
                c2.NextOnFailure = c3; c2.NextOnSuccess = c6;
                c3.NextOnFailure = c4; c3.NextOnSuccess = c5;
                c4.NextOnFailure = c4.NextOnSuccess = c5;
            
                // Once we find one on the second half-lattice, the rest are out.
                c5.NextOnFailure = c6; c5.NextOnSuccess = null;
                c6.NextOnFailure = c7; c6.NextOnSuccess = null;
                c7.NextOnFailure = c7.NextOnSuccess = null;
            
                LOOKUP_3D[i] = c0;
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
        }
    
        private class LatticePoint2D
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

        private class Grad2
        {
            public double dx, dy;
            public Grad2(double dx, double dy)
            {
                this.dx = dx; this.dy = dy;
            }
        }

        private class Grad3
        {
            public double dx, dy, dz;
            public Grad3(double dx, double dy, double dz)
            {
                this.dx = dx; this.dy = dy; this.dz = dz;
            }
        }
    }
}