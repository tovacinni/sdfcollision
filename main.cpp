#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/writeOBJ.h>
#include <igl/copyleft/marching_cubes.h>
#include <Eigen/Core>
#include "include/normalized_grid.h"
#include "include/normalized_cloth.h"
#include "include/distance.h"
#include "include/iqblob.h"

#include <iostream>

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
#define DEBUG_PRINT(x) std::cout << STRINGIFY(x) ":" << x << std::endl

int main(int argc, char *argv[])
{

    Eigen::MatrixXd Vcloth;
    Eigen::MatrixXi Fcloth;
    normalized_cloth(15, Vcloth, Fcloth);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;    

    const int res = 128;
    Eigen::MatrixXd GV;
    normalized_grid(res, GV);

    Eigen::VectorXd S(GV.rows());
    for (int i=0; i<GV.rows(); ++i) {
        S(i) = sdCone(GV.row(i));
    }

    igl::copyleft::marching_cubes(S,GV,res,res,res,V,F);
    
    Eigen::MatrixXd visG;
    normalized_grid(15, visG);

    Eigen::VectorXd visD(visG.rows());
    for (int i=0; i<visG.rows(); ++i) {
        visD(i) = sdCone(visG.row(i));
    }

    Eigen::MatrixXd inG;
    igl::slice_mask(visG, visD.array()<0.0, 1, inG);
    Eigen::MatrixXd outG;
    igl::slice_mask(visG, visD.array()>0.0, 1, outG);

    igl::opengl::glfw::Viewer viewer;
    std::cout<<R"(
  j       Move cloth downwards
  k       Move cloth upwards
)";

    viewer.data_list[0].set_mesh(V, F);
    viewer.data_list[0].set_points(outG, Eigen::RowVector3d(0.0,1.0,0.0));
    viewer.data_list[0].point_size = 3.0;
    viewer.data_list[0].show_faces = false;

    viewer.append_mesh();
    viewer.data_list[1].set_mesh(Vcloth,Fcloth);
    viewer.data_list[1].set_points(inG, Eigen::RowVector3d(1.0,0.0,0.0));
    //viewer.data_list[1].set_colors(Ccloth);
    viewer.data_list[1].point_size = 6.0;
    viewer.data_list[1].show_lines = false;
    viewer.data_list[1].face_based = false;
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = true;
    
    const auto & update = [&]() {
        Eigen::MatrixXd Ccloth(Fcloth.rows(), 3);
        for (int i=0; i<Fcloth.rows(); ++i) {
            double v0d = sdCone(Vcloth.row(Fcloth(i, 0)));
            double v1d = sdCone(Vcloth.row(Fcloth(i, 1)));
            double v2d = sdCone(Vcloth.row(Fcloth(i, 2)));
            if (v0d < 0.0 || v1d < 0.0 || v2d < 0.0) {
                Ccloth.row(i) = Eigen::RowVector3d(0.8,0.0,0.0);
            } else {
                Ccloth.row(i) = Eigen::RowVector3d(0.0,0.8,0.0);
            }
        }
        viewer.data_list[1].set_colors(Ccloth);

        viewer.data_list[1].set_vertices(Vcloth);
        viewer.data_list[1].compute_normals();
    };

    viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int){
        
        switch(key)
        {
            case 'j':
                Vcloth.col(1).array() -= 0.025;
                break;
            case 'k':
                Vcloth.col(1).array() += 0.025;
                break;
            default:
                return false;
        }
        update();
        return true;
    };

    
    viewer.launch();
    return EXIT_SUCCESS;
}

