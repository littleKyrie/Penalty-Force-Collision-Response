#include <V_linear_tetrahedron.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>
#include <iostream>

void V_linear_tetrahedron(double &energy,
                          Eigen::Ref<const Eigen::VectorXd> q,
                          Eigen::Ref<const Eigen::MatrixXd> V,
                          Eigen::Ref<const Eigen::RowVectorXi> element,
                          double volume, double C, double D) {

    // energy density formula
    auto neohookean_linear_tet = [&](
            double &e,
            Eigen::Ref<const Eigen::VectorXd>q,
            Eigen::Ref<const Eigen::Vector3d> X) {

        // compute dphi/dX
        Eigen::Matrix43d dphidX;
        dphi_linear_tetrahedron_dX(dphidX, V, element, X);

        // compute vertex matrix in world space
        Eigen::Matrix34d vertexMatrix;
        for (size_t j = 0; j < 4; j++) {
            auto colIndex = static_cast<Eigen::Index>(j);

            vertexMatrix(0, colIndex) = q[element[colIndex]*3+0];
            vertexMatrix(1, colIndex) = q[element[colIndex]*3+1];
            vertexMatrix(2, colIndex) = q[element[colIndex]*3+2];
        }

        // compute F
        Eigen::Matrix3d F;
        F = vertexMatrix * dphidX;

        // F is constant on the linear tet
        // it has the same integrated value on any point of the tet
        psi_neo_hookean(e, F, C, D); // energy density function on F
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);

    // integrate the energy density on the entire tet
    energy = energy * volume; // linear tet and constant F on X or tet
}