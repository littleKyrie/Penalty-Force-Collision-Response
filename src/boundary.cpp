//
// Created by Kyrie Zhang on 2023/12/16.
//

#include <boundary.h>

Boundary::Boundary() {
    point.setZero();
    // default as the floor
    normal[0] = 0.;
    normal[1] = 1.;
    normal[2] = 0.;
};

Boundary::Boundary(const Eigen::Vector3d& pos, const Eigen::Vector3d& dirToSimSpace) {
    point = pos;
    normal = dirToSimSpace;
}

void Boundary::setHashKeys(const std::vector<size_t> &keyList) {
    hashKeys = keyList;
}

