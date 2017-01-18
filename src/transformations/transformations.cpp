#include <math.h>
#include <transformations.h>
#include <iostream>
Matrix3d opkToRotation(float omega, float phi, float kappa) {
  Matrix3d o;
  o << 1, 0, 0,
       0, cos(omega), sin(omega),
       0, -sin(omega), cos(omega);

  Matrix3d p;
  p << cos(phi), 0, -sin(phi),
       0, 1, 0,
       sin(phi), 0, cos(phi);

  Matrix3d k;
  k << cos(kappa), sin(kappa), 0,
       -sin(kappa), cos(kappa), 0,
       0, 0, 1;

  //Chain multiplication roation = k*p*o
  k *= p;
  k *= o;
  return k;

}
