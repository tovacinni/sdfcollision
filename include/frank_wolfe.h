#ifndef FRANK_WOLFE_H
#define FRANK_WOLFE_H
#include <Eigen/Core>
#include "grad.h"

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
#define DEBUG_PRINT(x) std::cout << STRINGIFY(x) ":" << x << std::endl

// Given a triangle mesh and an implicit function, computes the minimum for each triangle.
// Inputs:
//   V  Matrix of the vertices of the mesh
//   F  Matrix of the faces of the mesh
//   f  Implicit function, which takes in a point and spits out a scalar
// Outputs:
//   P  Matrix of the minimum points
template<typename Sdf>
void frank_wolfe(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Sdf f,
    Eigen::MatrixXd & P) {

    P.resize(F.rows(), 3);
    for (int i=0; i<F.rows(); ++i) {
        Eigen::Vector3d v0 = V.row(F(i,0));
        Eigen::Vector3d v1 = V.row(F(i,1));
        Eigen::Vector3d v2 = V.row(F(i,2));
        Eigen::Vector3d x = (v0+v1+v2)/3.0;
#       pragma unroll
        for (int j=0; j<32; ++j) {
            double alpha = 1.0 / (j + 2.0);
            Eigen::Vector3d g;
            grad(x, f, g);
            double e0 = v0.dot(g);
            double e1 = v1.dot(g);
            double e2 = v2.dot(g);

            Eigen::Vector3d s = (e0<e1 && e0<e2) ? v0 : ((e1<e0 && e1<e2) ? v1 : v2);

            x += alpha * (s - x);
        }
        P.row(i) = x;
    }
    
}

#endif

