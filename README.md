# SuperSimplex & Fast Simplex-Style Gradient Noise

Successors to OpenSimplex Noise, plus updated OpenSimplex.

The provided 3D function in **SuperSimplexNoise** is about as fast as optimized OpenSimplex, but has better uniformity.

The provided 3D function in **FastSimplexStyleNoise** is about as fast as common Simplex noise implementations, but uses a much different process.

Both 2D functions are faster than the average.

All functions also include speed-optimized whole-area generators, which operate by flood-fill queue on the noise lattice.

All functions use new gradient sets that are symmetric with the lattice, but don't cause neighboring vertex gradients to constructively interfere.

Includes 2D and 3D noise. 4D noise is coming!


## Renders

### SuperSimplexNoise, 2D

![SuperSimplexNoise, 2D](images/ssn2.png?raw=true)

### FastSimplexStyleNoise, 2D

![FastSimplexStyleNoise, 2D](images/fssn2.png?raw=true)

### Updated OpenSimplexNoise, 2D

![SuperSimplexNoise, 3D (Classic, 2D slice)](images/ssn3c.png?raw=true)


### SuperSimplexNoise, 3D (Classic, 2D slice)

![Updated OpenSimplexNoise, 2D](images/osn2.png?raw=true)

### FastSimplexStyleNoise, 3D (Classic, 2D slice)

![FastSimplexStyleNoise, 3D (Classic, 2D slice)](images/fssn3c.png?raw=true)

### Updated OpenSimplexNoise, 3D (Classic, 2D slice)

![Updated OpenSimplexNoise, 3D (Classic, 2D slice)](images/osn3c.png?raw=true)


### SuperSimplexNoise, 3D (PlaneFirst, 2D slice)

![SuperSimplexNoise, 3D (PlaneFirst, 2D slice)](images/ssn3pf.png?raw=true)

### FastSimplexStyleNoise, 3D (PlaneFirst, 2D slice)

![FastSimplexStyleNoise, 3D (PlaneFirst, 2D slice)](images/fssn3pf.png?raw=true)

### Updated OpenSimplexNoise, 3D (PlaneFirst, 2D slices)

![Updated OpenSimplexNoise, 3D (PlaneFirst, 2D slice at z=0.0)](images/osn3pfa.png?raw=true)

![Updated OpenSimplexNoise, 3D (PlaneFirst, 2D slice at z=0.5)](images/osn3pfb.png?raw=true)


### Updated OpenSimplexNoise, 4D (2D slice)

![Updated OpenSimplexNoise, 4D (PlaneFirst, 2D slice)](images/osn4.png?raw=true)


## Performance Metrics

### SuperSimplex vs OpenSimplex (2D)

![2D Metrics SuperSimplexNoise](images/metrics_ssn2.png?raw=true)

It would appear that the older 2D OpenSimplex comes out ahead of the refactored version. But DigitalShadow's 3D and especially 4D refactorings come out ahead. The results may have been different had I tested in the original C# language of the refactored version, rather than Java.

### FastSimplexStyleNoise vs others (2D)

![2D Metrics FastSimplexStyleNoise](images/metrics_fssn2.png?raw=true)

### SuperSimplex vs OpenSimplex (3D)

![3D Metrics SuperSimplexNoise](images/metrics_ssn3.png?raw=true)

### FastSimplexStyleNoise vs others (3D)

![3D Metrics FastSimplexStyleNoise](images/metrics_fssn3.png?raw=true)

### OpenSimplex legacy vs DigitalShadow's refactor

![4D Metrics OpenSimplexNoise](images/metrics_osn4.png?raw=true)


## Public Domain Dedication

This is free and unencumbered software released into the public domain. See UNLICENSE. To the best of my non-lawyer knowledge, no patent claims cover anything implemented here (not legal advice). Please use this software to make cool things, rather than to make patents! Enjoy.