#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>
#include <iostream>

void d2V_linear_tetrahedron_dq2(
        Eigen::Matrix1212d &H,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::RowVectorXi> element,
        double volume, double C, double D) {

    auto neohookean_linear_tet = [&](
           Eigen::Matrix1212d &d2V,
           Eigen::Ref<const Eigen::VectorXd>q,
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

        //Code to compute non-integrated hessian matrix goes here
        Eigen::Matrix99d dP;
        d2psi_neo_hookean_dF2(dP, F, C, D);
        d2V = B.transpose() * dP * B;
    };

    //integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);

    H *= volume;
    

    //DO NOT REMOVE THIS CODE: This code ensures that the hessian matrix is symmetric positive definite by projecting all
    //negative eigenvalues to small, positive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
