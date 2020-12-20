#ifndef NORMALIZED_GRID_H
#define NORMALIZED_GRID_H
#include <Eigen/Core>

// Returns a normalized [-1,1] grid of coordinates for a volume.
// Inputs:
//   res # of vertices along each side
// Outputs:
//   GV  res^3 by 3 matrix of XYZ coords
void normalized_grid(
    const int res,
    Eigen::MatrixXd & GV);

#endif 

