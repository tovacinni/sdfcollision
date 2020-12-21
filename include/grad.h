#ifndef GRAD_H
#define GRAD_H
#include <Eigen/Core>

// Computes the gradient of an implicit function using finite differences.
// Inputs:
//   x  3d vector for the point at which to take the gradient
//   f  Implicit function, which takes in a point and spits out a scalar
// Outputs:
//   g  3d vector for the gradient
template<typename Sdf>
void grad(
    const Eigen::Vector3d & x,
    Sdf f,
    Eigen::Vector3d & g) {
    
    const double eps = 0.0003;
    const Eigen::Vector3d epx(eps, 0.0, 0.0);
    const Eigen::Vector3d epy(0.0, eps, 0.0);
    const Eigen::Vector3d epz(eps, 0.0, eps);

    g = Eigen::Vector3d(f(x+epx)-f(x-epx),
                        f(x+epy)-f(x-epy),
                        f(x+epz)-f(x-epz)).normalized();

}

#endif

