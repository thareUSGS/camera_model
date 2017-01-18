#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include <Eigen/Dense>

using namespace Eigen;

Matrix3d opkToRotation(float omega, float phi, float kappa);

#endif
