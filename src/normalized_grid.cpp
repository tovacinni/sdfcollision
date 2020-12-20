#include "../include/normalized_grid.h"

void normalized_grid(
    const int res,
    Eigen::MatrixXd & GV) {
    
    const Eigen::RowVector3i GVres(res, res, res);
    GV.resize(GVres(0)*GVres(1)*GVres(2), 3);
    for (int i=0; i<GVres(0); ++i) {
    for (int j=0; j<GVres(1); ++j) {
    for (int k=0; k<GVres(2); ++k) {
        double x = ((double) i / (double) GVres(0)) * 2.0 - 1.0;
        double y = ((double) j / (double) GVres(1)) * 2.0 - 1.0;
        double z = ((double) k / (double) GVres(2)) * 2.0 - 1.0;
        GV.row(i+GVres(0)*(j+GVres(1)*k)) = Eigen::RowVector3d(x,y,z);
    }}}

}

