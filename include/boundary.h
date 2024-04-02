//
// Created by Kyrie Zhang on 2023/12/12.
//

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <EigenTypes.h>
#include <Eigen/Dense>
#include <oneapi/tbb.h>
#include <soft_body.h>

// represent the plane boundary of the simulation
class Boundary {
public:
    Eigen::Vector3d point; // a point on the plane
    Eigen::Vector3d normal; // point to the simulation region
    std::vector<size_t> hashKeys;

public:
    Boundary();

    Boundary(const Eigen::Vector3d& pos, const Eigen::Vector3d& dirToSimSpace);

    void setHashKeys(const std::vector<size_t>& keyList);
};

typedef std::shared_ptr<Boundary> BoundaryPtr;

#endif //BOUNDARY_H
