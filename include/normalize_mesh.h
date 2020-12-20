#ifndef NORMALIZE_MESH_H
#define NORMALIZE_MESH_H
#include <Eigen/Core>

// Normalizes the mesh into a unit sphere.
// Inputs:
//   V  N by 3 matrix of vertices
//   F  N by 3 matrix of face indices
// Outputs:
//   V  N by 3 matrix of vertices
//   F  N by 3 matrix of face indices
void normalize_mesh(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & nV);
#endif
