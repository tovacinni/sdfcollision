#ifndef NORMALIZED_CLOTH_H
#define NORMALIZED_CLOTH_H
#include <Eigen/Core>

// Returns a normalized [-1,1] grid of coordinates and indices for a flat square cloth.
// Inputs:
//   res # of vertices along each side
// Outputs:
//   V  Matrix of vertices
//   F  Matrix of face indices
void normalized_cloth(
    const int res,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F);

#endif 

