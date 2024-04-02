//
// Created by Kyrie Zhang on 2023/12/16.
//

#include <math_tools.h>

Eigen::Vector4d barycentricCoordinate(
        const Eigen::Vector3d& node0,
        const Eigen::Vector3d& node1,
        const Eigen::Vector3d& node2,
        const Eigen::Vector3d& node3,
        const Eigen::Vector3d& position) {
    Eigen::Vector4d beta;

    Eigen::Vector3d X10 = node1 - node0;
    Eigen::Vector3d X20 = node2 - node0;
    Eigen::Vector3d X30 = node3 - node0;
    Eigen::Vector3d X12 = node1 - node2;
    Eigen::Vector3d X23 = node2 - node3;
    Eigen::Vector3d X13 = node1 - node3;
    Eigen::Vector3d X32 = - X23;
    Eigen::Vector3d X02 = - X20;
    Eigen::Vector3d X01 = - X10;
    Eigen::Vector3d X31 = - X13;
    Eigen::Vector3d Xp1 = position - node1;
    Eigen::Vector3d Xp2 = position - node2;
    Eigen::Vector3d Xp3 = position - node3;
    Eigen::Vector3d Xp0 = position - node0;

    // (0, 1, 2, 3)
    double v = X30.dot(X10.cross(X20)) / 6;
    // (3, 2, 1, p)
    double v0 = Xp3.dot(X23.cross(X13)) / 6;
    // (2, 3, 0, p)
    double v1 = Xp2.dot(X32.cross(X02)) / 6;
    // (1, 0, 3, p)
    double v2 = Xp1.dot(X01.cross(X31)) / 6;
    // (0, 1, 2, p)
    double v3 = Xp0.dot(X10.cross(X20)) / 6;

    beta[0] = v0 / v;
    beta[1] = v1 / v;
    beta[2] = v2 / v;
    beta[3] = v3 / v;

    return beta;
}

Eigen::Vector3d barycentricCoordinate(
        const Eigen::Vector3d& x0,
        const Eigen::Vector3d& x1,
        const Eigen::Vector3d& x2,
        const Eigen::Vector3d& p) {
    Eigen::Vector3d beta;
    Eigen::Vector3d localNormal; // this normal is only used to compute barycentric coordinate, may be opposite with the outward direction
    Eigen::Vector3d X10 = x1 - x0;
    Eigen::Vector3d X20 = x2 - x0;
    localNormal = X10.cross(X20);
    double area = localNormal.norm() / 2;
    localNormal = localNormal.normalized();
    Eigen::Vector3d X0p = x0 - p;
    Eigen::Vector3d X1p = x1 - p;
    Eigen::Vector3d X2p = x2 - p;
    double area0 = localNormal.dot(X1p.cross(X2p)) / 2;
    double area1 = localNormal.dot(X2p.cross(X0p)) / 2;
    double area2 = localNormal.dot(X0p.cross(X1p)) / 2;

    beta[0] = area0 / area;
    beta[1] = area1 / area;
    beta[2] = area2 / area;

    return beta;
}
