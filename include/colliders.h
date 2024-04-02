//
// Created by Kyrie Zhang on 2023/12/7.
//

#ifndef COLLIDERS_H
#define COLLIDERS_H

#include <soft_body.h>
#include <spatial_hash.h>
#include <EigenTypes.h>
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <testTools.h>

class Colliders {
    // data
public:
    std::vector<SoftBody> objects;

    std::vector<int> pointers; // conserve the initial index of q for every body

    std::vector<Boundary> boundaryList;

    // constructors
public:
    Colliders();

    explicit Colliders(const std::vector<SoftBody>& bodies);

    explicit Colliders(const SoftBody& object);

    Colliders(const std::vector<SoftBody>& bodies, const std::vector<int>& ptrs, const std::vector<Boundary>& boundaries);

    // uniform interfaces
public:
    virtual void collisionCull(const SpatialHash& hashTable, const Eigen::VectorXd& q);

    virtual void collisionForce(Eigen::VectorXd &f, const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot);

    virtual void boundaryForce(Eigen::VectorXd& f, const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot);

    virtual double collisionEnergy(SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot);

    // generate contact point
    virtual void computeConstraintSet(SpatialHash& hashTable, const Eigen::VectorXd& q);

    virtual void initialize();

    virtual void collisionStiffness(Eigen::SparseMatrixd &K, const Eigen::VectorXd &q);

    virtual void collisionDamping(Eigen::SparseMatrixd &K, const Eigen::VectorXd &qdot);
};

typedef std::shared_ptr<Colliders> CollidersPtr;

#endif //COLLIDERS_H
