## OpenSimplex 2 - SuperSimplex & Fast Simplex-Style Gradient Noise

Successors to OpenSimplex Noise, plus updated OpenSimplex. Includes 2D and 3D noise. 4D noise is coming!

* The provided 3D function in **OpenSimplex2S** (OpenSimplex 2, smooth version / "SuperSimplex") is about as fast as optimized OpenSimplex, but has better uniformity.

* The provided 3D function in **OpenSimplex2F** (OpenSimplex 2, faster version / "Fast Simplex-Style Noise") is about as fast as common Simplex noise implementations, but uses a much different process.

* The 2D functions aren't intended to represent new developments in the same vein as the 3D functions. They are just the logical pairings. Both 2D functions are implemented using lookup tables, use lattice-symmetric gradient sets, and perform similar to or faster than the average.

* All functions are given new gradient sets that are symmetric with the lattice, but don't cause neighboring vertex gradients to constructively interfere.

* Updated Legacy OpenSimplex with new gradient sets is included in java/legacy and csharp/legacy.

The classes in [java/areagen](https://github.com/KdotJPG/New-Simplex-Style-Gradient-Noise/tree/master/java/areagen) offer speed-optimized whole-area generators, which operate by flood-fill queue on the noise lattice. (i.e. they don't use a "range")

* The only differences between the two versions's area generators are: **radius parameter** and **normalization constant**.
  * The radius can be straightforwardly reduced for faster noise, or increased for smoother noise. There are no geometric traversal steps to cause discontinuities, and no hardcoded point computations to limit performance increases, as would be the case with varying the radius in the evaluators.
  * The normalization constant is baked into the same gradient set as the evaluator. It can be recomputed using [Noise Normalizer](https://github.com/KdotJPG/NoiseNormalizer). If left as is, the noise will still function correctly, it will just have a different output range.

#### TODO:

* 4D noise
* More language ports
* Consolidate render tiles in readme into fewer images
* Move radius into unified constant
* Test 24-gon instead of 12-gon for 2D gradient set, and replace if results nicer.
* Pull some of the explanation from [the reddit post](https://www.reddit.com/r/VoxelGameDev/comments/ee94wg/supersimplex_the_better_opensimplex_new_gradient/) into this readme.

#### Maybe TODO:

* Create combined OpenSimplex2F and OpenSimplex2S in one file, reducing repetition for someone who needs both.
* Include octave summation, ridged noise, etc.
* Exponentially-distributed noise ([source of idea](http://jcgt.org/published/0004/02/01/))
* Simultaneous multi-instance evaluation
* Disc/Ball-output noise (Outputs 2D/3D/etc. vector for more directionally-uniform domain warping)
* Tileable 2D noise (slightly mis-skewed triangular grid which repeats properly over a desired rectangle)
* Tileable 3D noise (using an extension of the above)
* Tileable 3D noise (exact, using the "classic" lattice orientation)

#### Change Log
* Replaced 12-direction 2D gradient set with a 24-direction set, to reduce visible feature repetition in thresholded single-octave 2D noise. (Feb 10, 2020)
* Renamed filenames FastSimplexStyleNoise to OpenSimplex2F, and SuperSimplexNoise to OpenSimplex2S. (Feb 10, 2020)
* Moved legacy OpenSimplex into legacy directories. (Feb 10, 2020)
* Slightly reorganized description above, and added TODO/changelog. (Jan 23, 2020)
* Renamed / additionally named the noise "OpenSimplex (2.0)", separated into two versions/variants. (Jan 23, 2020)
  * SuperSimplex and FastSimplexStyleNoise are very similar to each other algorithmically, and are in the same spirit as the original OpenSimplex.
  * OpenSimplex is used in a lot of projects, and the naming might help facilitate adoption of the new noise in its place.
* Add C# ports of evaluators. (Jan 13, 2020)
* Create separate files which include area generators. (Jan 13, 2020)
* Renamed PlaneFirst evaluators to XYBeforeZ, and added XZBeforeY. (Jan 13, 2020)
* Fixed equals() method in AreaGenLatticePoint3D for area generators. (Dec 24, 2019)

## Renders

### OpenSimplex2S, 2D

![OpenSimplex2S, 2D](images/ssn2.png?raw=true)

### OpenSimplex2F, 2D

![OpenSimplex2F, 2D](images/fssn2.png?raw=true)

### Updated OpenSimplexNoise, 2D

![OpenSimplex2S, 3D (Classic, 2D slice)](images/osn2.png?raw=true)

---

### OpenSimplex2S, 3D (Classic, 2D slice)

![Updated OpenSimplexNoise, 2D](images/ssn3c.png?raw=true)

### OpenSimplex2F, 3D (Classic, 2D slice)

![OpenSimplex2F, 3D (Classic, 2D slice)](images/fssn3c.png?raw=true)

### Updated OpenSimplexNoise, 3D (Classic, 2D slice)

![Updated OpenSimplexNoise, 3D (Classic, 2D slice)](images/osn3c.png?raw=true)

---

### OpenSimplex2S, 3D (XYBeforeZ, 2D slice)

![OpenSimplex2S, 3D (XYBeforeZ, 2D slice)](images/ssn3pf.png?raw=true)

### OpenSimplex2F, 3D (XYBeforeZ, 2D slice)

![OpenSimplex2F, 3D (XYBeforeZ, 2D slice)](images/fssn3pf.png?raw=true)

### Updated OpenSimplexNoise, 3D (XYBeforeZ, 2D slices)

![Updated OpenSimplexNoise, 3D (XYBeforeZ, 2D slice at z=0.0)](images/osn3pfa.png?raw=true)

![Updated OpenSimplexNoise, 3D (XYBeforeZ, 2D slice at z=0.5)](images/osn3pfb.png?raw=true)

---

### Updated OpenSimplex, 4D (2D slice)

![Updated OpenSimplex, 4D (2D slice)](images/osn4.png?raw=true)


## Performance Metrics

### SuperSimplex vs OpenSimplex (2D)

![2D Metrics OpenSimplex2S](images/metrics_ssn2.png?raw=true)

It would appear that the older 2D OpenSimplex comes out ahead of the refactored version. But DigitalShadow's 3D and especially 4D refactorings come out ahead. The results may have been different had I tested in the original C# language of the refactored version, rather than Java.

### OpenSimplex2F vs others (2D)

![2D Metrics OpenSimplex2F](images/metrics_fssn2.png?raw=true)

### SuperSimplex vs OpenSimplex (3D)

![3D Metrics OpenSimplex2S](images/metrics_ssn3.png?raw=true)

### OpenSimplex2F vs others (3D)

![3D Metrics OpenSimplex2F](images/metrics_fssn3.png?raw=true)

### OpenSimplex legacy vs DigitalShadow's refactor (4D)

![4D Metrics OpenSimplexNoise](images/metrics_osn4.png?raw=true)


## Public Domain Dedication

This is free and unencumbered software released into the public domain. See UNLICENSE. To the best of my non-lawyer knowledge, no patent claims cover anything implemented here (not legal advice). Please use this software to make cool things, rather than to make patents! Enjoy.
