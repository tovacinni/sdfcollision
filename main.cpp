#include <iostream>
#include <vector>

#include <Eigen/Core>

#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/writeOBJ.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/random_points_on_mesh.h>

#include "include/distance.h"
#include "include/frank_wolfe.h"
#include "include/grad.h"
#include "include/iqblob.h"
#include "include/normalized_cloth.h"
#include "include/normalized_grid.h"

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
#define DEBUG_PRINT(x) std::cout << STRINGIFY(x) ":" << x << std::endl

int main(int argc, char *argv[])
{
    // Build cloth
    Eigen::MatrixXd Vcloth;
    Eigen::MatrixXi Fcloth;
    normalized_cloth(15, Vcloth, Fcloth);
    Vcloth.col(0).array() += 0.05;
    //normalized_cloth(3, Vcloth, Fcloth);

    // Sample from cloth for sampling approach
    Eigen::MatrixXd Pcloth;
    Eigen::MatrixXi Icloth;
    Eigen::SparseMatrix<double> Bcloth;
    igl::random_points_on_mesh(5000, Vcloth, Fcloth, Bcloth, Icloth);
    Pcloth = Bcloth * Vcloth;

    Eigen::MatrixXd Pcloth_opt;
    frank_wolfe(Vcloth, Fcloth, sdCone, Pcloth_opt);

    for (int i=0; i<Vcloth.rows(); ++i) {
        Eigen::Vector3d g;
        grad(Vcloth.row(i), sdCone, g);
    }

    // Build Cone mesh with Marching Cubes
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;    

    const int res = 64;
    Eigen::MatrixXd GV;
    normalized_grid(res, GV);

    Eigen::VectorXd S(GV.rows());
    for (int i=0; i<GV.rows(); ++i) {
        S(i) = sdCone(GV.row(i));
    }

    igl::copyleft::marching_cubes(S,GV,res,res,res,V,F);
    
    // Build grid for visualization
    Eigen::MatrixXd visG;
    normalized_grid(8, visG);

    Eigen::VectorXd visD(visG.rows());
    for (int i=0; i<visG.rows(); ++i) {
        visD(i) = sdCone(visG.row(i));
    }

    Eigen::MatrixXd inG;
    igl::slice_mask(visG, visD.array()<0.0, 1, inG);
    Eigen::MatrixXd outG;
    igl::slice_mask(visG, visD.array()>0.0, 1, outG);

    // Building the visualization
    igl::opengl::glfw::Viewer viewer;
    std::cout<<R"(
  j       Move cloth downwards
  k       Move cloth upwards
)";
    
    viewer.core().background_color = Eigen::Vector4f(1.0, 1.0, 1.0, 1.0);

    viewer.data_list[0].set_mesh(V, F);
    viewer.data_list[0].point_size = 4.0;
    viewer.data_list[0].show_faces = false;
    viewer.data_list[0].line_width = 2.0;

    viewer.append_mesh();
    viewer.data_list[1].set_mesh(Vcloth,Fcloth);
    //viewer.data_list[1].set_colors(Ccloth);
    viewer.data_list[1].point_size = 8.0;
    viewer.data_list[1].line_width = 2.0;
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
            double od = sdCone(Pcloth_opt.row(i));
            if (od < 0.0 || v0d < 0.0 || v1d < 0.0 || v2d < 0.0) {
                Ccloth.row(i) = Eigen::RowVector3d(1.0,0.8,0.8);
            } else {
                Ccloth.row(i) = Eigen::RowVector3d(0.8,1.0,0.8);
            }
        }

        Eigen::VectorXd Vd(Vcloth.rows());
        for (int i=0; i<Vcloth.rows(); ++i) {
            Vd(i) = sdCone(Vcloth.row(i));
        }

        Eigen::MatrixXd inV;
        igl::slice_mask(Vcloth, Vd.array()<0.0, 1, inV);

        Eigen::VectorXd Pd(Pcloth.rows());
        for (int i=0; i<Pcloth.rows(); ++i) {
            Pd(i) = sdCone(Pcloth.row(i));
        }

        Eigen::MatrixXd inP;
        igl::slice_mask(Pcloth, Pd.array()<0.0, 1, inP);
        
        Eigen::MatrixXd outP;
        igl::slice_mask(Pcloth, Pd.array()>0.0, 1, outP);
        viewer.data_list[0].set_points(outG, Eigen::RowVector3d(0.0,1.0,0.0));
        viewer.data_list[0].add_points(inG, Eigen::RowVector3d(1.0,0.0,0.0));
        //viewer.data_list[0].add_points(inP, Eigen::RowVector3d(1.0,0.0,0.0));
        //viewer.data_list[0].add_points(outP, Eigen::RowVector3d(0.0,1.0,0.0));

        viewer.data_list[1].set_points(inV, Eigen::RowVector3d(1.0,0.0,0.0));
        frank_wolfe(Vcloth, Fcloth, sdCone, Pcloth_opt);
        viewer.data_list[1].add_points(Pcloth_opt, Eigen::RowVector3d(1.0,0.3,0.0));
        viewer.data_list[1].set_colors(Ccloth);

        viewer.data_list[1].set_vertices(Vcloth);
        viewer.data_list[1].compute_normals();
    };

    viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int){
        
        switch(key)
        {
            case 'j':
                //Vcloth.col(1).array() -= 0.025;
                Vcloth.col(1).array() -= 0.01;
                Pcloth.col(1).array() -= 0.01;
                break;
            case 'k':
                //Vcloth.col(1).array() += 0.025;
                Vcloth.col(1).array() += 0.01;
                Pcloth.col(1).array() += 0.01;
                break;
                break;
            default:
                return false;
        }
        update();
        return true;
    };

    update();
    viewer.launch();
    return EXIT_SUCCESS;
}

