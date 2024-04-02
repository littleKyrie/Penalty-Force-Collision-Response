#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(
        Eigen::Vector12d &dV,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::RowVectorXi> element,
        double volume, double C, double D) {
    /*
    // compute deformation gradient F
    // edge matrix T in reference space
    Eigen::Matrix3d T;
    Eigen::VectorXd X0 = V.row(element[0]);
    for (size_t i = 0; i < T.rows(); i++)
        for (size_t j = 0; j < T.cols(); j++) {
            // the col index of T will represent the response vertex index of the element
            // the row index of T will represent the response component of vertex coordinate
            auto vertexIndex = static_cast<Eigen::Index>(j);
            auto coordinate = static_cast<Eigen::Index>(i);

            T(coordinate, vertexIndex) =
                    V(element[vertexIndex+1], coordinate) - X0[coordinate];
        }
    */

    // in this function, it has the same structure as the V_linear_tetrahedron.cpp
    // but this part has nothing to do with quadrature
    auto neohookean_linear_tet = [&](
           Eigen::Vector12d &dV,
           Eigen::Ref<const Eigen::VectorXd>q,
           // Eigen::Ref<const Eigen::RowVectorXi> element,
           Eigen::Ref<const Eigen::Vector3d> X) {
        // dphi/dX
        Eigen::Matrix43d dphidX;
        dphi_linear_tetrahedron_dX(dphidX, V, element, X);

        // vertex matrix in world space
        Eigen::Matrix34d vertexMatrix;
        for (size_t j = 0; j < 4; j++) {
            auto colIndex = static_cast<Eigen::Index>(j);

            vertexMatrix(0, colIndex) = q[element[colIndex]*3+0];
            vertexMatrix(1, colIndex) = q[element[colIndex]*3+1];
            vertexMatrix(2, colIndex) = q[element[colIndex]*3+2];
        }

        // get F
        Eigen::Matrix3d F;
        F = vertexMatrix * dphidX;

        // compute the matrix B to transfer F from matrix to vector(row first)
        Eigen::MatrixXd B(9, 12);
        B.setZero();
        for (size_t j = 0; j < dphidX.rows(); j++) {
            auto colIndex = static_cast<Eigen::Index>(j);

            B(0, colIndex*3+0) = dphidX(colIndex, 0);
            B(1, colIndex*3+0) = dphidX(colIndex, 1);
            B(2, colIndex*3+0) = dphidX(colIndex, 2);
            B(3, colIndex*3+1) = dphidX(colIndex, 0);
            B(4, colIndex*3+1) = dphidX(colIndex, 1);
            B(5, colIndex*3+1) = dphidX(colIndex, 2);
            B(6, colIndex*3+2) = dphidX(colIndex, 0);
            B(7, colIndex*3+2) = dphidX(colIndex, 1);
            B(8, colIndex*3+2) = dphidX(colIndex, 2);
        }

        Eigen::Vector9d dpsi;
        dpsi_neo_hookean_dF(dpsi, F, C, D);
        dV = B.transpose() * dpsi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  

    dV *= volume;
}