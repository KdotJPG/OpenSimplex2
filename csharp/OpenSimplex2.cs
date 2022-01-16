/**
 * K.jpg's OpenSimplex 2, faster variant
 */

using System.Runtime.CompilerServices;

public static class OpenSimplex2
{
    private const long PRIME_X = 0x5205402B9270C86FL;
    private const long PRIME_Y = 0x598CD327003817B5L;
    private const long PRIME_Z = 0x5BCC226E9FA0BACBL;
    private const long PRIME_W = 0x56CC5227E58F554BL;
    private const long HASH_MULTIPLIER = 0x53A3F72DEEC546F5L;
    private const long SEED_FLIP_3D = -0x52D547B2E96ED629L;
    private const long SEED_OFFSET_4D = 0xE83DC3E0DA7164DL;

    private const double ROOT2OVER2 = 0.7071067811865476;
    private const double SKEW_2D = 0.366025403784439;
    private const double UNSKEW_2D = -0.21132486540518713;

    private const double ROOT3OVER3 = 0.577350269189626;
    private const double FALLBACK_ROTATE_3D = 2.0 / 3.0;
    private const double ROTATE_3D_ORTHOGONALIZER = UNSKEW_2D;

    private const float SKEW_4D = -0.138196601125011f;
    private const float UNSKEW_4D = 0.309016994374947f;
    private const float LATTICE_STEP_4D = 0.2f;

    private const int N_GRADS_2D_EXPONENT = 7;
    private const int N_GRADS_3D_EXPONENT = 8;
    private const int N_GRADS_4D_EXPONENT = 9;
    private const int N_GRADS_2D = 1 << N_GRADS_2D_EXPONENT;
    private const int N_GRADS_3D = 1 << N_GRADS_3D_EXPONENT;
    private const int N_GRADS_4D = 1 << N_GRADS_4D_EXPONENT;

    private const double NORMALIZER_2D = 0.01001634121365712;
    private const double NORMALIZER_3D = 0.07969837668935331;
    private const double NORMALIZER_4D = 0.0220065933241897;

    private const float RSQUARED_2D = 0.5f;
    private const float RSQUARED_3D = 0.6f;
    private const float RSQUARED_4D = 0.6f;


    /*
     * Noise Evaluators
     */

    /**
     * 2D Simplex noise, standard lattice orientation.
     */
    public static float Noise2(long seed, double x, double y)
    {
        // Get points for A2* lattice
        double s = SKEW_2D * (x + y);
        double xs = x + s, ys = y + s;

        return Noise2_UnskewedBase(seed, xs, ys);
    }

    /**
     * 2D Simplex noise, with Y pointing down the main diagonal.
     * Might be better for a 2D sandbox style game, where Y is vertical.
     * Probably slightly less optimal for heightmaps or continent maps,
     * unless your map is centered around an equator. It's a subtle
     * difference, but the option is here to make it an easy choice.
     */
    public static float Noise2_ImproveX(long seed, double x, double y)
    {
        // Skew transform and rotation baked into one.
        double xx = x * ROOT2OVER2;
        double yy = y * (ROOT2OVER2 * (1 + 2 * SKEW_2D));

        return Noise2_UnskewedBase(seed, yy + xx, yy - xx);
    }

    /**
     * 2D Simplex noise base.
     */
    private static float Noise2_UnskewedBase(long seed, double xs, double ys)
    {
        // Get base points and offsets.
        int xsb = FastFloor(xs), ysb = FastFloor(ys);
        float xi = (float)(xs - xsb), yi = (float)(ys - ysb);

        // Prime pre-multiplication for hash.
        long xsbp = xsb * PRIME_X, ysbp = ysb * PRIME_Y;

        // Unskew.
        float t = (xi + yi) * (float)UNSKEW_2D;
        float dx0 = xi + t, dy0 = yi + t;

        // First vertex.
        float value = 0;
        float a0 = RSQUARED_2D - dx0 * dx0 - dy0 * dy0;
        if (a0 > 0)
        {
            value = (a0 * a0) * (a0 * a0) * Grad(seed, xsbp, ysbp, dx0, dy0);
        }

        // Second vertex.
        float a1 = (float)(2 * (1 + 2 * UNSKEW_2D) * (1 / UNSKEW_2D + 2)) * t + ((float)(-2 * (1 + 2 * UNSKEW_2D) * (1 + 2 * UNSKEW_2D)) + a0);
        if (a1 > 0)
        {
            float dx1 = dx0 - (float)(1 + 2 * UNSKEW_2D);
            float dy1 = dy0 - (float)(1 + 2 * UNSKEW_2D);
            value += (a1 * a1) * (a1 * a1) * Grad(seed, xsbp + PRIME_X, ysbp + PRIME_Y, dx1, dy1);
        }

        // Third vertex.
        if (dy0 > dx0)
        {
            float dx2 = dx0 - (float)UNSKEW_2D;
            float dy2 = dy0 - (float)(UNSKEW_2D + 1);
            float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if (a2 > 0)
            {
                value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp, ysbp + PRIME_Y, dx2, dy2);
            }
        }
        else
        {
            float dx2 = dx0 - (float)(UNSKEW_2D + 1);
            float dy2 = dy0 - (float)UNSKEW_2D;
            float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if (a2 > 0)
            {
                value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp + PRIME_X, ysbp, dx2, dy2);
            }
        }

        return value;
    }

    /**
     * 3D OpenSimplex2 noise, with better visual isotropy in (X, Y).
     * Recommended for 3D terrain and time-varied animations.
     * The Z coordinate should always be the "different" coordinate in whatever your use case is.
     * If Y is vertical in world coordinates, call Noise3_ImproveXZ(x, z, Y) or use noise3_XZBeforeY.
     * If Z is vertical in world coordinates, call Noise3_ImproveXZ(x, y, Z).
     * For a time varied animation, call Noise3_ImproveXY(x, y, T).
     */
    public static float Noise3_ImproveXY(long seed, double x, double y, double z)
    {
        // Re-orient the cubic lattices without skewing, so Z points up the main lattice diagonal,
        // and the planes formed by XY are moved far out of alignment with the cube faces.
        // Orthonormal rotation. Not a skew transform.
        double xy = x + y;
        double s2 = xy * ROTATE_3D_ORTHOGONALIZER;
        double zz = z * ROOT3OVER3;
        double xr = x + s2 + zz;
        double yr = y + s2 + zz;
        double zr = xy * -ROOT3OVER3 + zz;

        // Evaluate both lattices to form a BCC lattice.
        return Noise3_UnrotatedBase(seed, xr, yr, zr);
    }

    /**
     * 3D OpenSimplex2 noise, with better visual isotropy in (X, Z).
     * Recommended for 3D terrain and time-varied animations.
     * The Y coordinate should always be the "different" coordinate in whatever your use case is.
     * If Y is vertical in world coordinates, call Noise3_ImproveXZ(x, Y, z).
     * If Z is vertical in world coordinates, call Noise3_ImproveXZ(x, Z, y) or use Noise3_ImproveXY.
     * For a time varied animation, call Noise3_ImproveXZ(x, T, y) or use Noise3_ImproveXY.
     */
    public static float Noise3_ImproveXZ(long seed, double x, double y, double z)
    {
        // Re-orient the cubic lattices without skewing, so Y points up the main lattice diagonal,
        // and the planes formed by XZ are moved far out of alignment with the cube faces.
        // Orthonormal rotation. Not a skew transform.
        double xz = x + z;
        double s2 = xz * ROTATE_3D_ORTHOGONALIZER;
        double yy = y * ROOT3OVER3;
        double xr = x + s2 + yy;
        double zr = z + s2 + yy;
        double yr = xz * -ROOT3OVER3 + yy;

        // Evaluate both lattices to form a BCC lattice.
        return Noise3_UnrotatedBase(seed, xr, yr, zr);
    }

    /**
     * 3D OpenSimplex2 noise, fallback rotation option
     * Use Noise3_ImproveXY or Noise3_ImproveXZ instead, wherever appropriate.
     * They have less diagonal bias. This function's best use is as a fallback.
     */
    public static float Noise3_Fallback(long seed, double x, double y, double z)
    {
        // Re-orient the cubic lattices via rotation, to produce a familiar look.
        // Orthonormal rotation. Not a skew transform.
        double r = FALLBACK_ROTATE_3D * (x + y + z);
        double xr = r - x, yr = r - y, zr = r - z;

        // Evaluate both lattices to form a BCC lattice.
        return Noise3_UnrotatedBase(seed, xr, yr, zr);
    }

    /**
     * Generate overlapping cubic lattices for 3D OpenSimplex2 noise.
     */
    private static float Noise3_UnrotatedBase(long seed, double xr, double yr, double zr)
    {
        // Get base points and offsets.
        int xrb = FastRound(xr), yrb = FastRound(yr), zrb = FastRound(zr);
        float xri = (float)(xr - xrb), yri = (float)(yr - yrb), zri = (float)(zr - zrb);

        // -1 if positive, 1 if negative.
        int xNSign = (int)(-1.0f - xri) | 1, yNSign = (int)(-1.0f - yri) | 1, zNSign = (int)(-1.0f - zri) | 1;

        // Compute absolute values, using the above as a shortcut. This was faster in my tests for some reason.
        float ax0 = xNSign * -xri, ay0 = yNSign * -yri, az0 = zNSign * -zri;

        // Prime pre-multiplication for hash.
        long xrbp = xrb * PRIME_X, yrbp = yrb * PRIME_Y, zrbp = zrb * PRIME_Z;

        // Loop: Pick an edge on each lattice copy.
        float value = 0;
        float a = (RSQUARED_3D - xri * xri) - (yri * yri + zri * zri);
        for (int l = 0; ; l++)
        {

            // Closest point on cube.
            if (a > 0)
            {
                value += (a * a) * (a * a) * Grad(seed, xrbp, yrbp, zrbp, xri, yri, zri);
            }

            // Second-closest point.
            if (ax0 >= ay0 && ax0 >= az0)
            {
                float b = a + ax0 + ax0;
                if (b > 1)
                {
                    b -= 1;
                    value += (b * b) * (b * b) * Grad(seed, xrbp - xNSign * PRIME_X, yrbp, zrbp, xri + xNSign, yri, zri);
                }
            }
            else if (ay0 > ax0 && ay0 >= az0)
            {
                float b = a + ay0 + ay0;
                if (b > 1)
                {
                    b -= 1;
                    value += (b * b) * (b * b) * Grad(seed, xrbp, yrbp - yNSign * PRIME_Y, zrbp, xri, yri + yNSign, zri);
                }
            }
            else
            {
                float b = a + az0 + az0;
                if (b > 1)
                {
                    b -= 1;
                    value += (b * b) * (b * b) * Grad(seed, xrbp, yrbp, zrbp - zNSign * PRIME_Z, xri, yri, zri + zNSign);
                }
            }

            // Break from loop if we're done, skipping updates below.
            if (l == 1) break;

            // Update absolute value.
            ax0 = 0.5f - ax0;
            ay0 = 0.5f - ay0;
            az0 = 0.5f - az0;

            // Update relative coordinate.
            xri = xNSign * ax0;
            yri = yNSign * ay0;
            zri = zNSign * az0;

            // Update falloff.
            a += (0.75f - ax0) - (ay0 + az0);

            // Update prime for hash.
            xrbp += (xNSign >> 1) & PRIME_X;
            yrbp += (yNSign >> 1) & PRIME_Y;
            zrbp += (zNSign >> 1) & PRIME_Z;

            // Update the reverse sign indicators.
            xNSign = -xNSign;
            yNSign = -yNSign;
            zNSign = -zNSign;

            // And finally update the seed for the other lattice copy.
            seed ^= SEED_FLIP_3D;
        }

        return value;
    }

    /**
     * 4D OpenSimplex2 noise, with XYZ oriented like Noise3_ImproveXY
     * and W for an extra degree of freedom. W repeats eventually.
     * Recommended for time-varied animations which texture a 3D object (W=time)
     * in a space where Z is vertical
     */
    public static float Noise4_ImproveXYZ_ImproveXY(long seed, double x, double y, double z, double w)
    {
        double xy = x + y;
        double s2 = xy * -0.21132486540518699998;
        double zz = z * 0.28867513459481294226;
        double ww = w * 0.2236067977499788;
        double xr = x + (zz + ww + s2), yr = y + (zz + ww + s2);
        double zr = xy * -0.57735026918962599998 + (zz + ww);
        double wr = z * -0.866025403784439 + ww;

        return Noise4_UnskewedBase(seed, xr, yr, zr, wr);
    }

    /**
     * 4D OpenSimplex2 noise, with XYZ oriented like Noise3_ImproveXZ
     * and W for an extra degree of freedom. W repeats eventually.
     * Recommended for time-varied animations which texture a 3D object (W=time)
     * in a space where Y is vertical
     */
    public static float Noise4_ImproveXYZ_ImproveXZ(long seed, double x, double y, double z, double w)
    {
        double xz = x + z;
        double s2 = xz * -0.21132486540518699998;
        double yy = y * 0.28867513459481294226;
        double ww = w * 0.2236067977499788;
        double xr = x + (yy + ww + s2), zr = z + (yy + ww + s2);
        double yr = xz * -0.57735026918962599998 + (yy + ww);
        double wr = y * -0.866025403784439 + ww;

        return Noise4_UnskewedBase(seed, xr, yr, zr, wr);
    }

    /**
     * 4D OpenSimplex2 noise, with XYZ oriented like Noise3_Fallback
     * and W for an extra degree of freedom. W repeats eventually.
     * Recommended for time-varied animations which texture a 3D object (W=time)
     * where there isn't a clear distinction between horizontal and vertical
     */
    public static float Noise4_ImproveXYZ(long seed, double x, double y, double z, double w)
    {
        double xyz = x + y + z;
        double ww = w * 0.2236067977499788;
        double s2 = xyz * -0.16666666666666666 + ww;
        double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

        return Noise4_UnskewedBase(seed, xs, ys, zs, ws);
    }

    /**
     * 4D OpenSimplex2 noise, fallback lattice orientation.
     */
    public static float Noise4_Fallback(long seed, double x, double y, double z, double w)
    {
        // Get points for A4 lattice
        double s = SKEW_4D * (x + y + z + w);
        double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

        return Noise4_UnskewedBase(seed, xs, ys, zs, ws);
    }

    /**
    * 4D OpenSimplex2 noise base.
    */
    private static float Noise4_UnskewedBase(long seed, double xs, double ys, double zs, double ws)
    {
        // Get base points and offsets
        int xsb = FastFloor(xs), ysb = FastFloor(ys), zsb = FastFloor(zs), wsb = FastFloor(ws);
        float xsi = (float)(xs - xsb), ysi = (float)(ys - ysb), zsi = (float)(zs - zsb), wsi = (float)(ws - wsb);

        // Determine which lattice we can be confident has a contributing point its corresponding cell's base simplex.
        // We only look at the spaces between the diagonal planes. This proved effective in all of my tests.
        float siSum = (xsi + ysi) + (zsi + wsi);
        int startingLattice = (int)(siSum * 1.25);

        // Offset for seed based on first lattice copy.
        seed += startingLattice * SEED_OFFSET_4D;

        // Offset for lattice point relative positions (skewed)
        float startingLatticeOffset = startingLattice * -LATTICE_STEP_4D;
        xsi += startingLatticeOffset; ysi += startingLatticeOffset; zsi += startingLatticeOffset; wsi += startingLatticeOffset;

        // Prep for vertex contributions.
        float ssi = (siSum + startingLatticeOffset * 4) * UNSKEW_4D;
        
        // Prime pre-multiplication for hash.
        long xsvp = xsb * PRIME_X, ysvp = ysb * PRIME_Y, zsvp = zsb * PRIME_Z, wsvp = wsb * PRIME_W;

        // Five points to add, total, from five copies of the A4 lattice.
        float value = 0;
        for (int i = 0; ; i++)
        {

            // Next point is the closest vertex on the 4-simplex whose base vertex is the aforementioned vertex.
            double score0 = 1.0 + ssi * (-1.0 / UNSKEW_4D); // Seems slightly faster than 1.0-xsi-ysi-zsi-wsi
            if (xsi >= ysi && xsi >= zsi && xsi >= wsi && xsi >= score0)
            {
                xsvp += PRIME_X;
                xsi -= 1;
                ssi -= UNSKEW_4D;
            }
            else if (ysi > xsi && ysi >= zsi && ysi >= wsi && ysi >= score0)
            {
                ysvp += PRIME_Y;
                ysi -= 1;
                ssi -= UNSKEW_4D;
            }
            else if (zsi > xsi && zsi > ysi && zsi >= wsi && zsi >= score0)
            {
                zsvp += PRIME_Z;
                zsi -= 1;
                ssi -= UNSKEW_4D;
            }
            else if (wsi > xsi && wsi > ysi && wsi > zsi && wsi >= score0)
            {
                wsvp += PRIME_W;
                wsi -= 1;
                ssi -= UNSKEW_4D;
            }

            // Gradient contribution with falloff.
            float dx = xsi + ssi, dy = ysi + ssi, dz = zsi + ssi, dw = wsi + ssi;
            float a = (dx * dx + dy * dy) + (dz * dz + dw * dw);
            if (a < RSQUARED_4D)
            {
                a -= RSQUARED_4D;
                a *= a;
                value += a * a * Grad(seed, xsvp, ysvp, zsvp, wsvp, dx, dy, dz, dw);
            }

            // Break from loop if we're done, skipping updates below.
            if (i == 4) break;

            // Update for next lattice copy shifted down by <-0.2, -0.2, -0.2, -0.2>.
            xsi += LATTICE_STEP_4D; ysi += LATTICE_STEP_4D; zsi += LATTICE_STEP_4D; wsi += LATTICE_STEP_4D;
            ssi += LATTICE_STEP_4D * 4 * UNSKEW_4D;
            seed -= SEED_OFFSET_4D;

            // Because we don't always start on the same lattice copy, there's a special reset case.
            if (i == startingLattice)
            {
                xsvp -= PRIME_X;
                ysvp -= PRIME_Y;
                zsvp -= PRIME_Z;
                wsvp -= PRIME_W;
                seed += SEED_OFFSET_4D * 5;
            }
        }

        return value;
    }

    /*
     * Utility
     */

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static float Grad(long seed, long xsvp, long ysvp, float dx, float dy)
    {
        long hash = seed ^ xsvp ^ ysvp;
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - N_GRADS_2D_EXPONENT + 1);
        int gi = (int)hash & ((N_GRADS_2D - 1) << 1);
        return GRADIENTS_2D[gi | 0] * dx + GRADIENTS_2D[gi | 1] * dy;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static float Grad(long seed, long xrvp, long yrvp, long zrvp, float dx, float dy, float dz)
    {
        long hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - N_GRADS_3D_EXPONENT + 2);
        int gi = (int)hash & ((N_GRADS_3D - 1) << 2);
        return GRADIENTS_3D[gi | 0] * dx + GRADIENTS_3D[gi | 1] * dy + GRADIENTS_3D[gi | 2] * dz;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static float Grad(long seed, long xsvp, long ysvp, long zsvp, long wsvp, float dx, float dy, float dz, float dw)
    {
        long hash = seed ^ (xsvp ^ ysvp) ^ (zsvp ^ wsvp);
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - N_GRADS_4D_EXPONENT + 2);
        int gi = (int)hash & ((N_GRADS_4D - 1) << 2);
        return (GRADIENTS_4D[gi | 0] * dx + GRADIENTS_4D[gi | 1] * dy) + (GRADIENTS_4D[gi | 2] * dz + GRADIENTS_4D[gi | 3] * dw);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int FastFloor(double x)
    {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int FastRound(double x)
    {
        return x < 0 ? (int)(x - 0.5) : (int)(x + 0.5);
    }

    /*
     * Gradients
     */

    private static readonly float[] GRADIENTS_2D;
    private static readonly float[] GRADIENTS_3D;
    private static readonly float[] GRADIENTS_4D;
    static OpenSimplex2()
    {

        GRADIENTS_2D = new float[N_GRADS_2D * 2];
        float[] grad2 = {
             0.38268343236509f,   0.923879532511287f,
             0.923879532511287f,  0.38268343236509f,
             0.923879532511287f, -0.38268343236509f,
             0.38268343236509f,  -0.923879532511287f,
            -0.38268343236509f,  -0.923879532511287f,
            -0.923879532511287f, -0.38268343236509f,
            -0.923879532511287f,  0.38268343236509f,
            -0.38268343236509f,   0.923879532511287f,
            //-------------------------------------//
             0.130526192220052f,  0.99144486137381f,
             0.608761429008721f,  0.793353340291235f,
             0.793353340291235f,  0.608761429008721f,
             0.99144486137381f,   0.130526192220051f,
             0.99144486137381f,  -0.130526192220051f,
             0.793353340291235f, -0.60876142900872f,
             0.608761429008721f, -0.793353340291235f,
             0.130526192220052f, -0.99144486137381f,
            -0.130526192220052f, -0.99144486137381f,
            -0.608761429008721f, -0.793353340291235f,
            -0.793353340291235f, -0.608761429008721f,
            -0.99144486137381f,  -0.130526192220052f,
            -0.99144486137381f,   0.130526192220051f,
            -0.793353340291235f,  0.608761429008721f,
            -0.608761429008721f,  0.793353340291235f,
            -0.130526192220052f,  0.99144486137381f,
        };
        for (int i = 0; i < grad2.Length; i++)
        {
            grad2[i] = (float)(grad2[i] / NORMALIZER_2D);
        }
        for (int i = 0, j = 0; i < GRADIENTS_2D.Length; i++, j++)
        {
            if (j == grad2.Length) j = 0;
            GRADIENTS_2D[i] = grad2[j];
        }

        GRADIENTS_3D = new float[N_GRADS_3D * 4];
        float[] grad3 = {
             2.22474487139f,       2.22474487139f,      -1.0f,                 0.0f,
             2.22474487139f,       2.22474487139f,       1.0f,                 0.0f,
             3.0862664687972017f,  1.1721513422464978f,  0.0f,                 0.0f,
             1.1721513422464978f,  3.0862664687972017f,  0.0f,                 0.0f,
            -2.22474487139f,       2.22474487139f,      -1.0f,                 0.0f,
            -2.22474487139f,       2.22474487139f,       1.0f,                 0.0f,
            -1.1721513422464978f,  3.0862664687972017f,  0.0f,                 0.0f,
            -3.0862664687972017f,  1.1721513422464978f,  0.0f,                 0.0f,
            -1.0f,                -2.22474487139f,      -2.22474487139f,       0.0f,
             1.0f,                -2.22474487139f,      -2.22474487139f,       0.0f,
             0.0f,                -3.0862664687972017f, -1.1721513422464978f,  0.0f,
             0.0f,                -1.1721513422464978f, -3.0862664687972017f,  0.0f,
            -1.0f,                -2.22474487139f,       2.22474487139f,       0.0f,
             1.0f,                -2.22474487139f,       2.22474487139f,       0.0f,
             0.0f,                -1.1721513422464978f,  3.0862664687972017f,  0.0f,
             0.0f,                -3.0862664687972017f,  1.1721513422464978f,  0.0f,
            //--------------------------------------------------------------------//
            -2.22474487139f,      -2.22474487139f,      -1.0f,                 0.0f,
            -2.22474487139f,      -2.22474487139f,       1.0f,                 0.0f,
            -3.0862664687972017f, -1.1721513422464978f,  0.0f,                 0.0f,
            -1.1721513422464978f, -3.0862664687972017f,  0.0f,                 0.0f,
            -2.22474487139f,      -1.0f,                -2.22474487139f,       0.0f,
            -2.22474487139f,       1.0f,                -2.22474487139f,       0.0f,
            -1.1721513422464978f,  0.0f,                -3.0862664687972017f,  0.0f,
            -3.0862664687972017f,  0.0f,                -1.1721513422464978f,  0.0f,
            -2.22474487139f,      -1.0f,                 2.22474487139f,       0.0f,
            -2.22474487139f,       1.0f,                 2.22474487139f,       0.0f,
            -3.0862664687972017f,  0.0f,                 1.1721513422464978f,  0.0f,
            -1.1721513422464978f,  0.0f,                 3.0862664687972017f,  0.0f,
            -1.0f,                 2.22474487139f,      -2.22474487139f,       0.0f,
             1.0f,                 2.22474487139f,      -2.22474487139f,       0.0f,
             0.0f,                 1.1721513422464978f, -3.0862664687972017f,  0.0f,
             0.0f,                 3.0862664687972017f, -1.1721513422464978f,  0.0f,
            -1.0f,                 2.22474487139f,       2.22474487139f,       0.0f,
             1.0f,                 2.22474487139f,       2.22474487139f,       0.0f,
             0.0f,                 3.0862664687972017f,  1.1721513422464978f,  0.0f,
             0.0f,                 1.1721513422464978f,  3.0862664687972017f,  0.0f,
             2.22474487139f,      -2.22474487139f,      -1.0f,                 0.0f,
             2.22474487139f,      -2.22474487139f,       1.0f,                 0.0f,
             1.1721513422464978f, -3.0862664687972017f,  0.0f,                 0.0f,
             3.0862664687972017f, -1.1721513422464978f,  0.0f,                 0.0f,
             2.22474487139f,      -1.0f,                -2.22474487139f,       0.0f,
             2.22474487139f,       1.0f,                -2.22474487139f,       0.0f,
             3.0862664687972017f,  0.0f,                -1.1721513422464978f,  0.0f,
             1.1721513422464978f,  0.0f,                -3.0862664687972017f,  0.0f,
             2.22474487139f,      -1.0f,                 2.22474487139f,       0.0f,
             2.22474487139f,       1.0f,                 2.22474487139f,       0.0f,
             1.1721513422464978f,  0.0f,                 3.0862664687972017f,  0.0f,
             3.0862664687972017f,  0.0f,                 1.1721513422464978f,  0.0f,
        };
        for (int i = 0; i < grad3.Length; i++)
        {
            grad3[i] = (float)(grad3[i] / NORMALIZER_3D);
        }
        for (int i = 0, j = 0; i < GRADIENTS_3D.Length; i++, j++)
        {
            if (j == grad3.Length) j = 0;
            GRADIENTS_3D[i] = grad3[j];
        }

        GRADIENTS_4D = new float[N_GRADS_4D * 4];
        float[] grad4 = {
            -0.6740059517812944f,   -0.3239847771997537f,   -0.3239847771997537f,    0.5794684678643381f,
            -0.7504883828755602f,   -0.4004672082940195f,    0.15296486218853164f,   0.5029860367700724f,
            -0.7504883828755602f,    0.15296486218853164f,  -0.4004672082940195f,    0.5029860367700724f,
            -0.8828161875373585f,    0.08164729285680945f,   0.08164729285680945f,   0.4553054119602712f,
            -0.4553054119602712f,   -0.08164729285680945f,  -0.08164729285680945f,   0.8828161875373585f,
            -0.5029860367700724f,   -0.15296486218853164f,   0.4004672082940195f,    0.7504883828755602f,
            -0.5029860367700724f,    0.4004672082940195f,   -0.15296486218853164f,   0.7504883828755602f,
            -0.5794684678643381f,    0.3239847771997537f,    0.3239847771997537f,    0.6740059517812944f,
            -0.6740059517812944f,   -0.3239847771997537f,    0.5794684678643381f,   -0.3239847771997537f,
            -0.7504883828755602f,   -0.4004672082940195f,    0.5029860367700724f,    0.15296486218853164f,
            -0.7504883828755602f,    0.15296486218853164f,   0.5029860367700724f,   -0.4004672082940195f,
            -0.8828161875373585f,    0.08164729285680945f,   0.4553054119602712f,    0.08164729285680945f,
            -0.4553054119602712f,   -0.08164729285680945f,   0.8828161875373585f,   -0.08164729285680945f,
            -0.5029860367700724f,   -0.15296486218853164f,   0.7504883828755602f,    0.4004672082940195f,
            -0.5029860367700724f,    0.4004672082940195f,    0.7504883828755602f,   -0.15296486218853164f,
            -0.5794684678643381f,    0.3239847771997537f,    0.6740059517812944f,    0.3239847771997537f,
            -0.6740059517812944f,    0.5794684678643381f,   -0.3239847771997537f,   -0.3239847771997537f,
            -0.7504883828755602f,    0.5029860367700724f,   -0.4004672082940195f,    0.15296486218853164f,
            -0.7504883828755602f,    0.5029860367700724f,    0.15296486218853164f,  -0.4004672082940195f,
            -0.8828161875373585f,    0.4553054119602712f,    0.08164729285680945f,   0.08164729285680945f,
            -0.4553054119602712f,    0.8828161875373585f,   -0.08164729285680945f,  -0.08164729285680945f,
            -0.5029860367700724f,    0.7504883828755602f,   -0.15296486218853164f,   0.4004672082940195f,
            -0.5029860367700724f,    0.7504883828755602f,    0.4004672082940195f,   -0.15296486218853164f,
            -0.5794684678643381f,    0.6740059517812944f,    0.3239847771997537f,    0.3239847771997537f,
             0.5794684678643381f,   -0.6740059517812944f,   -0.3239847771997537f,   -0.3239847771997537f,
             0.5029860367700724f,   -0.7504883828755602f,   -0.4004672082940195f,    0.15296486218853164f,
             0.5029860367700724f,   -0.7504883828755602f,    0.15296486218853164f,  -0.4004672082940195f,
             0.4553054119602712f,   -0.8828161875373585f,    0.08164729285680945f,   0.08164729285680945f,
             0.8828161875373585f,   -0.4553054119602712f,   -0.08164729285680945f,  -0.08164729285680945f,
             0.7504883828755602f,   -0.5029860367700724f,   -0.15296486218853164f,   0.4004672082940195f,
             0.7504883828755602f,   -0.5029860367700724f,    0.4004672082940195f,   -0.15296486218853164f,
             0.6740059517812944f,   -0.5794684678643381f,    0.3239847771997537f,    0.3239847771997537f,
            //------------------------------------------------------------------------------------------//
            -0.753341017856078f,    -0.37968289875261624f,  -0.37968289875261624f,  -0.37968289875261624f,
            -0.7821684431180708f,   -0.4321472685365301f,   -0.4321472685365301f,    0.12128480194602098f,
            -0.7821684431180708f,   -0.4321472685365301f,    0.12128480194602098f,  -0.4321472685365301f,
            -0.7821684431180708f,    0.12128480194602098f,  -0.4321472685365301f,   -0.4321472685365301f,
            -0.8586508742123365f,   -0.508629699630796f,     0.044802370851755174f,  0.044802370851755174f,
            -0.8586508742123365f,    0.044802370851755174f, -0.508629699630796f,     0.044802370851755174f,
            -0.8586508742123365f,    0.044802370851755174f,  0.044802370851755174f, -0.508629699630796f,
            -0.9982828964265062f,   -0.03381941603233842f,  -0.03381941603233842f,  -0.03381941603233842f,
            -0.37968289875261624f,  -0.753341017856078f,    -0.37968289875261624f,  -0.37968289875261624f,
            -0.4321472685365301f,   -0.7821684431180708f,   -0.4321472685365301f,    0.12128480194602098f,
            -0.4321472685365301f,   -0.7821684431180708f,    0.12128480194602098f,  -0.4321472685365301f,
             0.12128480194602098f,  -0.7821684431180708f,   -0.4321472685365301f,   -0.4321472685365301f,
            -0.508629699630796f,    -0.8586508742123365f,    0.044802370851755174f,  0.044802370851755174f,
             0.044802370851755174f, -0.8586508742123365f,   -0.508629699630796f,     0.044802370851755174f,
             0.044802370851755174f, -0.8586508742123365f,    0.044802370851755174f, -0.508629699630796f,
            -0.03381941603233842f,  -0.9982828964265062f,   -0.03381941603233842f,  -0.03381941603233842f,
            -0.37968289875261624f,  -0.37968289875261624f,  -0.753341017856078f,    -0.37968289875261624f,
            -0.4321472685365301f,   -0.4321472685365301f,   -0.7821684431180708f,    0.12128480194602098f,
            -0.4321472685365301f,    0.12128480194602098f,  -0.7821684431180708f,   -0.4321472685365301f,
             0.12128480194602098f,  -0.4321472685365301f,   -0.7821684431180708f,   -0.4321472685365301f,
            -0.508629699630796f,     0.044802370851755174f, -0.8586508742123365f,    0.044802370851755174f,
             0.044802370851755174f, -0.508629699630796f,    -0.8586508742123365f,    0.044802370851755174f,
             0.044802370851755174f,  0.044802370851755174f, -0.8586508742123365f,   -0.508629699630796f,
            -0.03381941603233842f,  -0.03381941603233842f,  -0.9982828964265062f,   -0.03381941603233842f,
            -0.37968289875261624f,  -0.37968289875261624f,  -0.37968289875261624f,  -0.753341017856078f,
            -0.4321472685365301f,   -0.4321472685365301f,    0.12128480194602098f,  -0.7821684431180708f,
            -0.4321472685365301f,    0.12128480194602098f,  -0.4321472685365301f,   -0.7821684431180708f,
             0.12128480194602098f,  -0.4321472685365301f,   -0.4321472685365301f,   -0.7821684431180708f,
            -0.508629699630796f,     0.044802370851755174f,  0.044802370851755174f, -0.8586508742123365f,
             0.044802370851755174f, -0.508629699630796f,     0.044802370851755174f, -0.8586508742123365f,
             0.044802370851755174f,  0.044802370851755174f, -0.508629699630796f,    -0.8586508742123365f,
            -0.03381941603233842f,  -0.03381941603233842f,  -0.03381941603233842f,  -0.9982828964265062f,
            -0.3239847771997537f,   -0.6740059517812944f,   -0.3239847771997537f,    0.5794684678643381f,
            -0.4004672082940195f,   -0.7504883828755602f,    0.15296486218853164f,   0.5029860367700724f,
             0.15296486218853164f,  -0.7504883828755602f,   -0.4004672082940195f,    0.5029860367700724f,
             0.08164729285680945f,  -0.8828161875373585f,    0.08164729285680945f,   0.4553054119602712f,
            -0.08164729285680945f,  -0.4553054119602712f,   -0.08164729285680945f,   0.8828161875373585f,
            -0.15296486218853164f,  -0.5029860367700724f,    0.4004672082940195f,    0.7504883828755602f,
             0.4004672082940195f,   -0.5029860367700724f,   -0.15296486218853164f,   0.7504883828755602f,
             0.3239847771997537f,   -0.5794684678643381f,    0.3239847771997537f,    0.6740059517812944f,
            -0.3239847771997537f,   -0.3239847771997537f,   -0.6740059517812944f,    0.5794684678643381f,
            -0.4004672082940195f,    0.15296486218853164f,  -0.7504883828755602f,    0.5029860367700724f,
             0.15296486218853164f,  -0.4004672082940195f,   -0.7504883828755602f,    0.5029860367700724f,
             0.08164729285680945f,   0.08164729285680945f,  -0.8828161875373585f,    0.4553054119602712f,
            -0.08164729285680945f,  -0.08164729285680945f,  -0.4553054119602712f,    0.8828161875373585f,
            -0.15296486218853164f,   0.4004672082940195f,   -0.5029860367700724f,    0.7504883828755602f,
             0.4004672082940195f,   -0.15296486218853164f,  -0.5029860367700724f,    0.7504883828755602f,
             0.3239847771997537f,    0.3239847771997537f,   -0.5794684678643381f,    0.6740059517812944f,
            -0.3239847771997537f,   -0.6740059517812944f,    0.5794684678643381f,   -0.3239847771997537f,
            -0.4004672082940195f,   -0.7504883828755602f,    0.5029860367700724f,    0.15296486218853164f,
             0.15296486218853164f,  -0.7504883828755602f,    0.5029860367700724f,   -0.4004672082940195f,
             0.08164729285680945f,  -0.8828161875373585f,    0.4553054119602712f,    0.08164729285680945f,
            -0.08164729285680945f,  -0.4553054119602712f,    0.8828161875373585f,   -0.08164729285680945f,
            -0.15296486218853164f,  -0.5029860367700724f,    0.7504883828755602f,    0.4004672082940195f,
             0.4004672082940195f,   -0.5029860367700724f,    0.7504883828755602f,   -0.15296486218853164f,
             0.3239847771997537f,   -0.5794684678643381f,    0.6740059517812944f,    0.3239847771997537f,
            -0.3239847771997537f,   -0.3239847771997537f,    0.5794684678643381f,   -0.6740059517812944f,
            -0.4004672082940195f,    0.15296486218853164f,   0.5029860367700724f,   -0.7504883828755602f,
             0.15296486218853164f,  -0.4004672082940195f,    0.5029860367700724f,   -0.7504883828755602f,
             0.08164729285680945f,   0.08164729285680945f,   0.4553054119602712f,   -0.8828161875373585f,
            -0.08164729285680945f,  -0.08164729285680945f,   0.8828161875373585f,   -0.4553054119602712f,
            -0.15296486218853164f,   0.4004672082940195f,    0.7504883828755602f,   -0.5029860367700724f,
             0.4004672082940195f,   -0.15296486218853164f,   0.7504883828755602f,   -0.5029860367700724f,
             0.3239847771997537f,    0.3239847771997537f,    0.6740059517812944f,   -0.5794684678643381f,
            -0.3239847771997537f,    0.5794684678643381f,   -0.6740059517812944f,   -0.3239847771997537f,
            -0.4004672082940195f,    0.5029860367700724f,   -0.7504883828755602f,    0.15296486218853164f,
             0.15296486218853164f,   0.5029860367700724f,   -0.7504883828755602f,   -0.4004672082940195f,
             0.08164729285680945f,   0.4553054119602712f,   -0.8828161875373585f,    0.08164729285680945f,
            -0.08164729285680945f,   0.8828161875373585f,   -0.4553054119602712f,   -0.08164729285680945f,
            -0.15296486218853164f,   0.7504883828755602f,   -0.5029860367700724f,    0.4004672082940195f,
             0.4004672082940195f,    0.7504883828755602f,   -0.5029860367700724f,   -0.15296486218853164f,
             0.3239847771997537f,    0.6740059517812944f,   -0.5794684678643381f,    0.3239847771997537f,
            -0.3239847771997537f,    0.5794684678643381f,   -0.3239847771997537f,   -0.6740059517812944f,
            -0.4004672082940195f,    0.5029860367700724f,    0.15296486218853164f,  -0.7504883828755602f,
             0.15296486218853164f,   0.5029860367700724f,   -0.4004672082940195f,   -0.7504883828755602f,
             0.08164729285680945f,   0.4553054119602712f,    0.08164729285680945f,  -0.8828161875373585f,
            -0.08164729285680945f,   0.8828161875373585f,   -0.08164729285680945f,  -0.4553054119602712f,
            -0.15296486218853164f,   0.7504883828755602f,    0.4004672082940195f,   -0.5029860367700724f,
             0.4004672082940195f,    0.7504883828755602f,   -0.15296486218853164f,  -0.5029860367700724f,
             0.3239847771997537f,    0.6740059517812944f,    0.3239847771997537f,   -0.5794684678643381f,
             0.5794684678643381f,   -0.3239847771997537f,   -0.6740059517812944f,   -0.3239847771997537f,
             0.5029860367700724f,   -0.4004672082940195f,   -0.7504883828755602f,    0.15296486218853164f,
             0.5029860367700724f,    0.15296486218853164f,  -0.7504883828755602f,   -0.4004672082940195f,
             0.4553054119602712f,    0.08164729285680945f,  -0.8828161875373585f,    0.08164729285680945f,
             0.8828161875373585f,   -0.08164729285680945f,  -0.4553054119602712f,   -0.08164729285680945f,
             0.7504883828755602f,   -0.15296486218853164f,  -0.5029860367700724f,    0.4004672082940195f,
             0.7504883828755602f,    0.4004672082940195f,   -0.5029860367700724f,   -0.15296486218853164f,
             0.6740059517812944f,    0.3239847771997537f,   -0.5794684678643381f,    0.3239847771997537f,
             0.5794684678643381f,   -0.3239847771997537f,   -0.3239847771997537f,   -0.6740059517812944f,
             0.5029860367700724f,   -0.4004672082940195f,    0.15296486218853164f,  -0.7504883828755602f,
             0.5029860367700724f,    0.15296486218853164f,  -0.4004672082940195f,   -0.7504883828755602f,
             0.4553054119602712f,    0.08164729285680945f,   0.08164729285680945f,  -0.8828161875373585f,
             0.8828161875373585f,   -0.08164729285680945f,  -0.08164729285680945f,  -0.4553054119602712f,
             0.7504883828755602f,   -0.15296486218853164f,   0.4004672082940195f,   -0.5029860367700724f,
             0.7504883828755602f,    0.4004672082940195f,   -0.15296486218853164f,  -0.5029860367700724f,
             0.6740059517812944f,    0.3239847771997537f,    0.3239847771997537f,   -0.5794684678643381f,
             0.03381941603233842f,   0.03381941603233842f,   0.03381941603233842f,   0.9982828964265062f,
            -0.044802370851755174f, -0.044802370851755174f,  0.508629699630796f,     0.8586508742123365f,
            -0.044802370851755174f,  0.508629699630796f,    -0.044802370851755174f,  0.8586508742123365f,
            -0.12128480194602098f,   0.4321472685365301f,    0.4321472685365301f,    0.7821684431180708f,
             0.508629699630796f,    -0.044802370851755174f, -0.044802370851755174f,  0.8586508742123365f,
             0.4321472685365301f,   -0.12128480194602098f,   0.4321472685365301f,    0.7821684431180708f,
             0.4321472685365301f,    0.4321472685365301f,   -0.12128480194602098f,   0.7821684431180708f,
             0.37968289875261624f,   0.37968289875261624f,   0.37968289875261624f,   0.753341017856078f,
             0.03381941603233842f,   0.03381941603233842f,   0.9982828964265062f,    0.03381941603233842f,
            -0.044802370851755174f,  0.044802370851755174f,  0.8586508742123365f,    0.508629699630796f,
            -0.044802370851755174f,  0.508629699630796f,     0.8586508742123365f,   -0.044802370851755174f,
            -0.12128480194602098f,   0.4321472685365301f,    0.7821684431180708f,    0.4321472685365301f,
             0.508629699630796f,    -0.044802370851755174f,  0.8586508742123365f,   -0.044802370851755174f,
             0.4321472685365301f,   -0.12128480194602098f,   0.7821684431180708f,    0.4321472685365301f,
             0.4321472685365301f,    0.4321472685365301f,    0.7821684431180708f,   -0.12128480194602098f,
             0.37968289875261624f,   0.37968289875261624f,   0.753341017856078f,     0.37968289875261624f,
             0.03381941603233842f,   0.9982828964265062f,    0.03381941603233842f,   0.03381941603233842f,
            -0.044802370851755174f,  0.8586508742123365f,   -0.044802370851755174f,  0.508629699630796f,
            -0.044802370851755174f,  0.8586508742123365f,    0.508629699630796f,    -0.044802370851755174f,
            -0.12128480194602098f,   0.7821684431180708f,    0.4321472685365301f,    0.4321472685365301f,
             0.508629699630796f,     0.8586508742123365f,   -0.044802370851755174f, -0.044802370851755174f,
             0.4321472685365301f,    0.7821684431180708f,   -0.12128480194602098f,   0.4321472685365301f,
             0.4321472685365301f,    0.7821684431180708f,    0.4321472685365301f,   -0.12128480194602098f,
             0.37968289875261624f,   0.753341017856078f,     0.37968289875261624f,   0.37968289875261624f,
             0.9982828964265062f,    0.03381941603233842f,   0.03381941603233842f,   0.03381941603233842f,
             0.8586508742123365f,   -0.044802370851755174f, -0.044802370851755174f,  0.508629699630796f,
             0.8586508742123365f,   -0.044802370851755174f,  0.508629699630796f,    -0.044802370851755174f,
             0.7821684431180708f,   -0.12128480194602098f,   0.4321472685365301f,    0.4321472685365301f,
             0.8586508742123365f,    0.508629699630796f,    -0.044802370851755174f, -0.044802370851755174f,
             0.7821684431180708f,    0.4321472685365301f,   -0.12128480194602098f,   0.4321472685365301f,
             0.7821684431180708f,    0.4321472685365301f,    0.4321472685365301f,   -0.12128480194602098f,
             0.753341017856078f,     0.37968289875261624f,   0.37968289875261624f,   0.37968289875261624f,
        };
        for (int i = 0; i < grad4.Length; i++)
        {
            grad4[i] = (float)(grad4[i] / NORMALIZER_4D);
        }
        for (int i = 0, j = 0; i < GRADIENTS_4D.Length; i++, j++)
        {
            if (j == grad4.Length) j = 0;
            GRADIENTS_4D[i] = grad4[j];
        }
    }
}
