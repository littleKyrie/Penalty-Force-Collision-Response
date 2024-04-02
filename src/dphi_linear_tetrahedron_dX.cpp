#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>

void dphi_linear_tetrahedron_dX(
        Eigen::Matrix43d &dphi,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::RowVectorXi> element,
        Eigen::Ref<const Eigen::Vector3d> X) {

    Eigen::Matrix3d T;
    Eigen::VectorXd X0 = V.row(element[0]);
    Eigen::VectorXd X1 = V.row(element[1]);
    Eigen::VectorXd X2 = V.row(element[2]);
    Eigen::VectorXd X3 = V.row(element[3]);

    T.block<3, 1>(0, 0) = X1 - X0;
    T.block<3, 1>(0, 1) = X2 - X0;
    T.block<3, 1>(0, 2) = X3 - X0;

    Eigen::RowVector3d one;
    one << 1.0, 1.0, 1.0;

    dphi.setZero();
    dphi.block<1, 3>(0, 0) = -one * T.inverse();
    dphi.block<3, 3>(1, 0) = T.inverse();
}