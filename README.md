# OpenSimplex 2

Successors to OpenSimplex Noise, plus updated OpenSimplex. Includes 2D, 3D, and 4D noise.

## Motivation

Legacy OpenSimplex used a different grid layout scheme to Simplex. It worked fine for common applications, but didn't produce as consistent contrast as the original Simplex point layout in 3D/4D. It also broke down completely when extended to 5D+, producing diagonal bands of visibly-varying value ranges.

Instead of vertex rearrangement, OpenSimplex2 focuses on congruent layouts constructed and traversed using disparate methods where need was concluded.

Gradient vector tables were also revisited to improve probability symmetry in both the old and new noise.

## Included

### OpenSimplex2(F)
 * Looks the most like Simplex.
 * Is about as fast as common Simplex implementations.
 * 3D is implemented by constructing a rotated body-centered-cubic grid as the offset union of two rotated cubic grids, and finding the closest two points on each.
 * 4D is implemented using 5 offset copies of the dual (reverse-skewed) grid to produce the target grid. This technique generalizes to any number of dimensions.
 * Bundled 2D function is just 2D Simplex with a nice gradient table. The 4D technique does extend to 2D, but this was found to be faster and I couldn't find evidence of limitations on its use. 3D and 4D are novel contributions to the field.

### OpenSimplex2S
 * Looks the most like 2014 OpenSimplex.
 * Uses large vertex contribution ranges like 2014 OpenSimplex, but has better uniformity in 3D and 4D.
 * 3D is implemented analogously to OpenSimplex2(F), but finds the closest four points on each to match the larger radius.
 * 4D is implemented using the ordinary skew and a pre-generated 4x4x4x4 lookup table. I have a work-in-progress more-procedural implementation in my project backlog, which I will swap in if I find it to be competitive in performance.
 * 2D is based on Simplex, but uses a different process to find which points are in range.
 * Recommended choice for ridged noise (if passing individual layers into `abs(x)`).

### Legacy OpenSimplex with updated gradients
 * In some code bases using original OpenSimplex, it may take fewer dev cycles to swap in the same noise that just has the new gradient tables. This brings moderate improvements without changing consistency characteristics or internal frequency.
 * Included in `_old/{language}/legacy`.

### Older OpenSimplex2 implementations
 * OpenSimplex2 noise implementations from before performance and portability re-engineering.
 * Can also be found in the `_old` directory.

Note: area-generators have been moved to [their original repository](https://github.com/KdotJPG/Noise-VertexQueue-AreaGen).

## Changelog

* Tuned up this `README.md`. (Mar 26, 2022)
* Re-wrote functions to be instancelessly seedable and less dependent on lookup tables. Re-organized repository. Renamed `OpenSimplex2F` to just `OpenSimplex2` in file/class names. (Jan 16, 2022)
* Shortened lookup table for Simplex/OpenSimplex2(F) 4D (July 5, 2020)
* Added 4D to OpenSimplex2(F) Java/C#, port OpenSimplex2S 4D to C#. (July 5, 2020)
* Added 4D to OpenSimplex2S Java. (Apr 30, 2020)
* Replaced individual renders in `README.md` with consolidated renders. (Feb 10, 2020)
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

* OpenSimplex2S is a smoother copy of OpenSimplex2(F).
* OpenSimplex2S and Updated OpenSimplex are effectively identical.
* Original OpenSimplex produced more straight parts and was not as probabilistically symmetric.

### 3D Noise (ImproveXY Orientation)

![Renders 3D ImproveXY](https://user-images.githubusercontent.com/8829856/149676980-768775b4-f0ca-417b-aa0a-2379b61fb69d.png)

* OpenSimplex2(F) and OpenSimplex2S, 2D slices of 3D ImproveXY, keep mostly the same properties as you move along the third axis.
* Updated OpenSimplex looks good in both slices, but the slices look different from each other. In an animation, this is particularly noticeable.
* Original OpenSimplex was less uniform and not as probabilistically symmetric.

### 3D Noise (Fallback Orientation)

![Renders 3D Fallback](https://user-images.githubusercontent.com/8829856/149676981-2d40fff9-b585-4fca-9a09-eb51d5927fec.png)

* OpenSimplex2(F) and OpenSimplex2S, 2D slices of 3D in the classic lattice orientation, look decent but are less ideal for X/Y planes being the primary focus.
* This Updated OpenSimplex render appears to show less directional bias than original OpenSimplex.

### 4D Noise (ImproveXY_ImproveZW Orientation)

![Renders 4D ImproveXY_ImproveZW](https://user-images.githubusercontent.com/8829856/149676999-14dc7dfc-a94c-4a8a-b654-00e25f33eb68.png)

* OpenSimplex2S 4D has higher apparent contrast than original OpenSimplex.
* OpenSimplex2(F) 4D has a very dotty appearance. Fine for fBm, but may pair best with extra steps if used for ridged noise.

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

This is free and unencumbered software. The majority of files are released under CC0. Where marked, files in certain directories fall under UNLICENSE instead, as they were based off of UNLICENSE code other than mine (e.g. Pull Requests, legacy OpenSimplex with DigitalShadow's lookup tables). To the best of my non-lawyer knowledge, no patent claims cover anything implemented here (though also nothing here is legal advice). Please use this software to make cool things, rather than to make patents! Enjoy.
