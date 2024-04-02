//
// Created by Kyrie Zhang on 2023/12/8.
//

#ifndef MATH_TOOLS_H
#define MATH_TOOLS_H

#include <EigenTypes.h>
#include <Eigen/Dense>

Eigen::Vector4d barycentricCoordinate(
        const Eigen::Vector3d& node0,
        const Eigen::Vector3d& node1,
        const Eigen::Vector3d& node2,
        const Eigen::Vector3d& node3,
        const Eigen::Vector3d& position);

Eigen::Vector3d barycentricCoordinate(
        const Eigen::Vector3d& x0,
        const Eigen::Vector3d& x1,
        const Eigen::Vector3d& x2,
        const Eigen::Vector3d& p);

#endif //MATH_TOOLS_H
