#ifndef GRAD_H
#define GRAD_H
#include <Eigen/Core>

typedef void (*callback_function)(void);

// Computes the gradient of an implicit function using finite differences.
// Inputs:
//   x  3d vector for the point at which to take the gradient
//   f  Implicit function, which takes in a point and spits out a scalar
// Outputs:
//   g  3d vector for the gradient
void grad(
    const Eigen::Vector3d & x,
    callback_function f,
    Eigen::Vector3d & g);

#endif

