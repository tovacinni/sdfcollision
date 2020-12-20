
#include "../include/normalized_grid.h"

void normalized_cloth(
    const int res,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F) {
    
    const double h = 1.0 / res;

    V.resize(res*res, 3);
    F.resize((res-1)*(res-1)*2, 3);

    for (int i=0; i<V.rows(); ++i) {
        int r = i / res;
        int c = i % res;
        V.row(i) = Eigen::RowVector3d(r*h, 0.5, c*h);
    }
    V.array() = V.array() * 2.0 - 1.0;

    for (int i=0; i<F.rows()/2; ++i) {
        int r = i / (res-1);
        int c = i % (res-1);
        F.row(2*i  ) = Eigen::RowVector3i(r*res+c, r*res+c+1, (r+1)*res+c);  
        F.row(2*i+1) = Eigen::RowVector3i((r+1)*res+c, r*res+c+1, (r+1)*res+c+1);  
    }

}

