#include "../include/grad.h"

void grad(
    const Eigen::Vector3d & x,
    callback_function f,
    Eigen::Vector3d & g) {
    
    const double eps = 0.0003;
    const Eigen::Vector3d epx(eps, 0.0, 0.0);
    const Eigen::Vector3d epy(0.0, eps, 0.0);
    const Eigen::Vector3d epz(eps, 0.0, eps);

    g = Eigen::Vector3d(f(x+epx)-f(x-epx),
                        f(x+epy)-f(x-epy),
                        f(x+epz)-f(x-epz)).normalized();

}

