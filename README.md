Adaptive Gaussian Quadrature with Modern Fortran

### Brief description

An object-oriented modern Fortran library to integrate functions using adaptive Gaussian quadrature. There are five selectable methods to use:
* Adaptive 6-point Legendre-Gauss
* Adaptive 8-point Legendre-Gauss
* Adaptive 10-point Legendre-Gauss
* Adaptive 12-point Legendre-Gauss
* Adaptive 14-point Legendre-Gauss

The core code is based on the SLATEC routine [DGAUS8](http://www.netlib.org/slatec/src/dgaus8.f) (which is the source of the 8-point routine). Coefficients for the others were obtained from [here](http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php).

### License

The quadrature-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/quadrature-module/blob/master/LICENSE) (BSD-style).

### Keywords

* adaptive quadrature, automatic integrator, gauss quadrature, numerical integration
