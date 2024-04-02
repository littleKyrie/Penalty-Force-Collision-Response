//
// Created by Kyrie Zhang on 2024/3/10.
//

#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <EigenTypes.h>
#include <Eigen/Dense>
#include <cmath>

// represent a simple rotation by matrix
void rotateAroundZ(Eigen::Matrix3d &R, double angle);

#endif //TRANSFORM_H
