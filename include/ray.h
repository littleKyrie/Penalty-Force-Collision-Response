//
// Created by Kyrie Zhang on 2023/12/16.
//

#ifndef RAY_H
#define RAY_H

#include <EigenTypes.h>
#include <Eigen/Dense>

// class of a single ray (more like an edge in my definition)
class Ray {
    // data
public:
    Eigen::Vector3d origin;
    Eigen::Vector3d dir; // ||dir|| = edge length
    double t; // impact of time

public:
    Ray();

    Ray(const Eigen::Vector3d& o, const Eigen::Vector3d& d);

public:
    // ray(edge) intersects with triangle by Moller-Trumbore algorithm
    bool intersectPoint(Eigen::Vector3d& w, const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);
};

#endif //RAY_H
