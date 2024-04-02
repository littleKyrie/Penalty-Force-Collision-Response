//
// Created by Kyrie Zhang on 2023/12/16.
//

#include <ray.h>

Ray::Ray() {
    origin.setZero();
    // default: x-axis positive direction
    dir.setZero();
    dir[0] = 1;
    t = std::numeric_limits<double>::max();
}

Ray::Ray(const Eigen::Vector3d& o, const Eigen::Vector3d& d) {
    origin = o;
    dir = d;
    t = 1.;
}

bool Ray::intersectPoint(
        Eigen::Vector3d& w,
        const Eigen::Vector3d& p0,
        const Eigen::Vector3d& p1,
        const Eigen::Vector3d& p2) {
    // p = alpha * p0 + beta * p1 + gamma * p2 (barycentric coordinate presentation)
    // solve equations: (-d, p1-p0, p2-p0)(t, beta, gamma)T = (o-p0)
    // t = ((o - p0) x (p1 - p0)) * (p2 - p0) / (d x (p2 - p0)) * (p1 - p0)
    // beta = (d x (p2 - p0)) * (o - p0) / (d x (p2 - p0)) * (p1 - p0)
    // gama = (d x (o - p0)) * (p1 - p0) / (d x (p2 - p0)) * (p1 - p0)
    Eigen::Vector3d b = origin - p0;
    Eigen::Vector3d p10 = p1 - p0;
    Eigen::Vector3d p20 = p2 - p0;
    double dno = p10.dot(dir.cross(p20));
    double r = p20.dot(b.cross(p10)) / dno;
    double beta = b.dot(dir.cross(p20)) / dno;
    double gamma = p10.dot(dir.cross(b)) / dno;
    double alpha = 1 - beta - gamma;

    if (r >= 0 && r < 1 && beta >= 0 && gamma >= 0 && alpha >= 0) {
        t = r;
        w[0] = alpha;
        w[1] = beta;
        w[2] = gamma;
        return true;
    }

    return false;
}
