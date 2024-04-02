#include <V_spring_particle_particle.h>

//the potential energy of a spring with 3D end points q0 and q1 and undeformed length l0
void V_spring_particle_particle(
        double &V,
        Eigen ::Ref<const Eigen::Vector3d> q0,
        Eigen::Ref<const Eigen::Vector3d> q1,
        double l0, double stiffness) {
    double lengthSquared = (q1 - q0).dot(q1 - q0);
    double length = sqrt(lengthSquared);
    V = 0.5 * stiffness * (length - l0) * (length - l0);
}