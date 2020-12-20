#ifndef DISTANCE_H
#define DISTANCE_H
#include <Eigen/Core>

// Distance functions taken from https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

double sdSphere(const Eigen::Vector3d & p) {
    return p.squaredNorm() - 0.5;
}

double sdCone(const Eigen::Vector3d & p) {
    double q = sqrt(pow(p(0), 2.0) + pow(p(2), 2.0));
    return std::max(sin(1.11)*q+cos(1.11)*p(1), -1.0-p(1));
}

#endif 

