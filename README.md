**quadrature-fortran** : Adaptive Gaussian Quadrature with Modern Fortran

### Brief description

An object-oriented modern Fortran library to integrate functions using adaptive Gaussian quadrature. There are five selectable methods to use:
* Adaptive 6-point Legendre-Gauss
* Adaptive 8-point Legendre-Gauss
* Adaptive 10-point Legendre-Gauss
* Adaptive 12-point Legendre-Gauss
* Adaptive 14-point Legendre-Gauss

The library supports:

* 1D integration:
  ```math
  \int_{x_l}^{x_u} f(x) dx
  ```
* 2D integration:
  ```math
  \int_{y_l}^{y_u} \int_{x_l}^{x_u} f(x,y) dx dy
  ```
* 3D integration:
  ```math
  \int_{z_l}^{z_u} \int_{y_l}^{y_u} \int_{x_l}^{x_u} f(x,y,z) dx dy dz
  ```
* 4D integration:
  ```math
  \int_{q_l}^{q_u} \int_{z_l}^{z_u} \int_{y_l}^{y_u} \int_{x_l}^{x_u} f(x,y,z,q) dx dy dz dq
  ```
* 5D integration:
  ```math
  \int_{r_l}^{r_u} \int_{q_l}^{q_u} \int_{z_l}^{z_u} \int_{y_l}^{y_u} \int_{x_l}^{x_u} f(x,y,z,q,r) dx dy dz dq dr
  ```
* 6D integration:
  ```math
  \int_{s_l}^{s_u} \int_{r_l}^{r_u} \int_{q_l}^{q_u} \int_{z_l}^{z_u} \int_{y_l}^{y_u} \int_{x_l}^{x_u} f(x,y,z,q,r,s) dx dy dz dq dr ds
  ```

The core code is based on the SLATEC routine [DGAUS8](http://www.netlib.org/slatec/src/dgaus8.f) (which is the source of the 8-point routine). Coefficients for the others were obtained from [here](http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php). The original 1D code has been generalized for multi-dimensional integration.

### License

The quadrature-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/quadrature-fortran/blob/master/LICENSE) (BSD-style).

### Keywords

* adaptive quadrature, automatic integrator, gauss quadrature, numerical integration
