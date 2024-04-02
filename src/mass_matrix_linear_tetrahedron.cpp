 
 #include <mass_matrix_linear_tetrahedron.h>

#include <iostream>

 void mass_matrix_linear_tetrahedron(
         Eigen::Matrix1212d &M,
         Eigen::Ref<const Eigen::VectorXd> qdot,
         Eigen::Ref<const Eigen::RowVectorXi> element,
         double density, double volume) {
    // compute mass matrix based on barycentric coordinates shape functions
    // compute the integral of the base functions production on the volume as the coefficient of mass matrix component
    double integralBase[4][4];
    integralBase[0][0] = 0.10f * density * volume;
    integralBase[0][1] = 0.05f * density * volume;
    integralBase[0][2] = 0.05f * density * volume;
    integralBase[0][3] = 0.05f * density * volume;
    integralBase[1][0] = 0.05f * density * volume;
    integralBase[1][1] = 0.10f * density * volume;
    integralBase[1][2] = 0.05f * density * volume;
    integralBase[1][3] = 0.05f * density * volume;
    integralBase[2][0] = 0.05f * density * volume;
    integralBase[2][1] = 0.05f * density * volume;
    integralBase[2][2] = 0.10f * density * volume;
    integralBase[2][3] = 0.05f * density * volume;
    integralBase[3][0] = 0.05f * density * volume;
    integralBase[3][1] = 0.05f * density * volume;
    integralBase[3][2] = 0.05f * density * volume;
    integralBase[3][3] = 0.10f * density * volume;

    M.setZero();
    // assemble the mass matrix
    for (size_t i = 0; i < 4; i++)
        for (size_t j = 0; j < 4; j++) {
            auto rowIndex = static_cast<Eigen::Index>(i);
            auto colIndex = static_cast<Eigen::Index>(j);

            M.block<3, 3>(3*rowIndex, 3*colIndex) = integralBase[i][j] * Eigen::Matrix3d::Identity();
        }

 }