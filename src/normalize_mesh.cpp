#include "../include/normalize_mesh.h"

#include <iostream>

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
#define DEBUG_PRINT(x) std::cout << STRINGIFY(x) ":" << x << std::endl

void normalize_mesh(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & nV) {
    
    Eigen::VectorXd max_V = V.colwise().maxCoeff();
    Eigen::VectorXd min_V = V.colwise().minCoeff();
    Eigen::RowVectorXd center = ((max_V + min_V) / 2.0);

    Eigen::MatrixXd centered_V = V.rowwise() - center;

    double max_norm = centered_V.rowwise().norm().maxCoeff();
    nV = centered_V * (1.0 / max_norm);
}

