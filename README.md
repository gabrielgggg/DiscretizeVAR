# DiscretizeVAR
Discretize VAR(1) of arbitrary size, with arbitrary covariance matrix for innovations. Support for VAR(1) with covariance matrix perturbed by common AR(1) volatility shock, e.g. "volatility regime," like baseline Bansal-Yaron. Allows the elimination of support points with low probability in the ergodic distribution (non-tensor grid).

Uses the Armadillo library for C++, with HDF5 support for I-O.

![Example plot](example_2d.png)
