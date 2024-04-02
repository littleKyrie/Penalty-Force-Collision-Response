//
// Created by Kyrie Zhang on 2023/12/7.
//

#include <colliders.h>

Colliders::Colliders() = default;

Colliders::Colliders(const SoftBody &object) {
    objects.push_back(object);
    pointers.push_back(0);
}

Colliders::Colliders(const std::vector<SoftBody> &colliders) {
    pointers.resize(colliders.size());
    pointers[0] = 0;
    objects.push_back(colliders[0]);
    for (size_t i = 1; i < colliders.size(); i++) {
        objects.push_back(colliders[i]);
        pointers[i] = pointers[i-1] + 3 * colliders[i-1].V.rows();
    }
}

Colliders::Colliders(
        const std::vector<SoftBody> &bodies,
        const std::vector<int> &ptrs,
        const std::vector<Boundary> &boundaries) {
    objects = bodies;
    pointers = ptrs;
    boundaryList = boundaries;
}

void Colliders::collisionCull(const SpatialHash &hashTable, const Eigen::VectorXd &q) {
}

void Colliders::collisionForce(Eigen::VectorXd &f, const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot) {
}

void Colliders::boundaryForce(Eigen::VectorXd &f, const SpatialHash &hashTable, const Eigen::VectorXd &q,
                              const Eigen::VectorXd &qdot) {
}

double Colliders::collisionEnergy(SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot) {
    return 0;
}

void Colliders::computeConstraintSet(SpatialHash &hashTable, const Eigen::VectorXd &q) {
}

void Colliders::initialize() {
}

void Colliders::collisionStiffness(Eigen::SparseMatrixd &K, const Eigen::VectorXd &q) {
}

void Colliders::collisionDamping(Eigen::SparseMatrixd &K, const Eigen::VectorXd &qdot) {
}
