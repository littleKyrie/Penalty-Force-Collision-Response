//
// Created by Kyrie Zhang on 2023/12/7.
//

#ifndef PENALTY_FORCE_COLLISION_H
#define PENALTY_FORCE_COLLISION_H

#include <colliders.h>
#include <ray.h>
#include <oneapi/tbb.h>

class PenaltyForceColliders : public Colliders {
    // typedef
public:
    struct IntersectionPoint {
        Eigen::Vector3d _normal; // the intersected surface unit outward normal
        double _t; // the distance scalar factor from the border point to the exact intersection point
        Eigen::Vector3d _depth; // the distance vector from the border point pointing to the exact intersection point

        IntersectionPoint();
    };

    struct ContactPoint {
        bool _colliding;
        bool _border;
        bool _processed;
        int _tet; // the collided tet number
        int _obj; // the collided tet belonged to which object
        Eigen::Vector4d _beta; // barycentric coordinates of a tetrahedron
        double _depth; // consistent depth to the collided surface
        Eigen::Vector3d _normal; // consistent direction of the spring force
        std::vector<size_t> _intersectNeighbors; // indices of adjacent intersecting edges (for intersectingEdgeSet or exactPoints)
        std::vector<size_t> _collidingNeighbors; // indices of adjacent colliding points (for collidingPoints)

        ContactPoint();

        void initialize();
    };

    // data
private:
    double _stiffness = 10;
    double _damping = 1e-3;
    std::vector<std::vector<ContactPoint>> _pointList; // contains all points
    std::vector<Eigen::MatrixXi> _intersectingEdgeSet; // extracted form the edge matrix for edge-surface intersection test
    std::vector<std::vector<IntersectionPoint>> _exactPoints; // stores intersecting points of every intersecting edges
    std::vector<std::vector<int>> _collidingPoints; // pointers point colliding points in the pointList;

    // constructors
public:
    PenaltyForceColliders();

    explicit PenaltyForceColliders(const std::vector<SoftBody>& colliders);

    explicit PenaltyForceColliders(const SoftBody &object);

    PenaltyForceColliders(const std::vector<SoftBody>& bodies, const std::vector<int>& ptrs, const std::vector<Boundary>& boundaries);

    // resolve collision
public:
    void collisionCull(const SpatialHash& hashTable, const Eigen::VectorXd& q) override;

    void collisionForce(Eigen::VectorXd &f, const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot) override;

    double collisionEnergy(SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot) override;

    void computeConstraintSet(SpatialHash &hashTable, const Eigen::VectorXd& q) override;

    void initialize() override;

    void collisionStiffness(Eigen::SparseMatrixd &K, const Eigen::VectorXd &q) override;

    void collisionDamping(Eigen::SparseMatrixd &K, const Eigen::VectorXd &qdot) override;

    void constructIntersectingEdgeSet();

    // this has been done in the collisionCull
    void constructCollidingPoints();

    void constructCollidingNeighbors();

    void computeExactIntersectionPoint(const SpatialHash& hashTable, const Eigen::VectorXd& q);

    void computeDepthOfBorderPoints();

    void propagate(const Eigen::VectorXd& q);

    double collisionResponse(Eigen::VectorXd &f, const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot);

    void assemblePenaltyForce(Eigen::VectorXd &f, const Eigen::VectorXd& qdot);

    double penaltyEnergy();

    double collideWithBoundary(Eigen::VectorXd& f, const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot);

    // divided from the function "collideWithBoundary"
    void boundaryForce(Eigen::VectorXd& f, const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot);

    // divided from the function "collideWithBoundary"
    double boundaryEnergy(const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot);

    // setters
public:
    void setStiffness(double k);

    void setDamping(double d);

    // for test
public:
    void visualColliders() const;
};

typedef std::shared_ptr<PenaltyForceColliders> PenaltyForceCollidersPtr;

#endif //PENALTY_FORCE_COLLISION_H