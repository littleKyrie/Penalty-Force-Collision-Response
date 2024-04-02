//
// Created by Kyrie Zhang on 2024/3/10.
//

#include <transform.h>

void rotateAroundZ(Eigen::Matrix3d &R, double angle) {
    R.setIdentity();

    double angleRadians = angle * M_PI / 180.0;
    double sine = std::sin(angleRadians);
    double cosine = std::cos(angleRadians);
    R.coeffRef(0, 0) = cosine;
    R.coeffRef(0, 1) = -sine;
    R.coeffRef(1, 0) = sine;
    R.coeffRef(1, 1) = cosine;
}