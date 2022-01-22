# OpenSimplex 2

Successors to OpenSimplex Noise, plus updated OpenSimplex. Includes 2D, 3D, and 4D noise.

* The provided 3D and 4D functions in **OpenSimplex2S** (smooth version / "SuperSimplex") are about as fast as optimized OpenSimplex, but have better uniformity.

* The provided 3D and 4D functions in **OpenSimplex2(F)** (faster version) are about as fast as common Simplex noise implementations, but use much different processes.

* The 2D function in OpenSimplex2(F) is just 2D Simplex with a nice gradient table. The rest are novel contributions to the field.

* All functions are given gradient sets that are symmetric with the lattice, but don't cause neighboring vertex gradients to constructively interfere.

* Updated Legacy OpenSimplex with new gradient sets is included in `_old/java/legacy` and `_old/csharp/legacy`. Older implementations of OpenSimplex2 can also be found in the `_old` directory.

Note: area-generators have been moved to [the original repository](https://github.com/KdotJPG/Noise-VertexQueue-AreaGen)

#### TODO:

* More language ports
* Pull some of the explanation from [the reddit post](https://www.reddit.com/r/VoxelGameDev/comments/ee94wg/supersimplex_the_better_opensimplex_new_gradient/) into this readme.

#### Maybe TODO:

* Create combined OpenSimplex2(F) and OpenSimplex2S in one file, reducing repetition for someone who needs both.
* Include octave summation (fBm), ridged noise, etc.
* Simultaneous multi-instance evaluation.
* Disc/Ball-output noise (Outputs 2D/3D/etc. vector for more directionally-uniform domain warping)
* Tileable 2D noise (slightly mis-skewed triangular grid which repeats properly over a desired rectangle)
* Tileable 3D noise (using an extension of the above)
* Tileable 3D noise (exact, using the "fallback" lattice orientation)

#### Change Log

* Re-wrote functions to be instancelessly seedable and less dependent on lookup tables. Re-organized repository. (Jan 16, 2021)
* Shortened lookup table for Simplex/OpenSimplex2F 4D (July 5, 2020)
* Add 4D to OpenSimplex2F Java/C#, port OpenSimplex2S 4D to C#. (July 5, 2020)
* Add 4D to OpenSimplex2S Java. (Apr 30, 2020)
* Replaced individual renders in readme with consolidated renders. (Feb 10, 2020)
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

### 2D Noise

![Renders 2D](https://user-images.githubusercontent.com/8829856/149676975-ba7cfc82-3bb5-4e0c-b6cb-8315eca51193.png)

* OpenSimplex2S is a smoother copy of OpenSimplex2F.
* OpenSimplex2S and Updated OpenSimplex are effectively identical.
* Original OpenSimplex produced more straight parts and was not probabilistically lattice-symmetric.

### 3D Noise (ImproveXY Orientation)

![Renders 3D ImproveXY](https://user-images.githubusercontent.com/8829856/149676980-768775b4-f0ca-417b-aa0a-2379b61fb69d.png)

* OpenSimplex2(F) and OpenSimplex2S, 2D slices of 3D ImproveXY, keep mostly the same properties as you move along the third axis.
* Updated OpenSimplex looks good in both slices, but the slices look different from each other. In an animation, this is particularly noticeable.
* Original OpenSimplex was less uniform and not probabilistically lattice-symmetric.

### 3D Noise (Fallback Orientation)

![Renders 3D Fallback](https://user-images.githubusercontent.com/8829856/149676981-2d40fff9-b585-4fca-9a09-eb51d5927fec.png)

* OpenSimplex2 and OpenSimplex2S, 2D slices of 3D in the classic lattice orientation, look decent but are less ideal for X/Y planes being the primary focus.
* This Updated OpenSimplex render appears to show less directional bias than original OpenSimplex.

### 4D Noise (ImproveXY_ImproveZW Orientation)

![Renders 4D ImproveXY_ImproveZW](https://user-images.githubusercontent.com/8829856/149676999-14dc7dfc-a94c-4a8a-b654-00e25f33eb68.png)

* OpenSimplex2S 4D has higher apparent contrast than original OpenSimplex.
* OpenSimplex2(F) 4D has a very dotty appearance. Fine for fBm but probably not ridged noise.

### 4D Noise (ImproveXYZ Orientation)

![Renders 4D ImproveXYZ](https://user-images.githubusercontent.com/8829856/149677002-542aefc5-3e60-46d2-ac4a-6d947c0b94a7.png)

* 2D slices look fine, but this rotation is most intended for texturing 3D objects with an additional time variable.

### 4D Noise (Fallback Orientation)

![Renders 4D Fallback](https://user-images.githubusercontent.com/8829856/149677011-267ab94c-3f3b-47bf-85a2-7b00befa95d6.png)

* OpenSimplex2(F) has the most noticeable diagonal artifacts, followed by Old OpenSimplex.

### 4D Noise (Torus Mapping)

![Renders 4D Torus](https://user-images.githubusercontent.com/8829856/149677013-c0593b51-b757-41ed-9ece-1bf2f88feee7.png)

* Seamless tileable 2D noise from 4D, mapped using noise4(r sin x, r cos x, r sin y, r cos y)

## Public Domain Dedication

This is free and unencumbered software. The majority of files are released under CC0. Where marked, files in certain directories fall under UNLICENSE instead, as they came from Pull Requests when that was the active license for the repository (or were otherwise based off of UNLICENSE code, i.e. legacy OpenSimplex with DigitalShadow's lookup tables). To the best of my non-lawyer knowledge, no patent claims cover anything implemented here (not legal advice). Please use this software to make cool things, rather than to make patents! Enjoy.
