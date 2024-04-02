//
// Created by Kyrie Zhang on 2023/12/7.
//

#include <penalty_force_collision.h>
#include <math_tools.h>

PenaltyForceColliders::ContactPoint::ContactPoint() {
    _colliding = false;
    _border = false;
    _processed = false;
    _tet = -1;
    _obj = -1;
    _beta.setZero();
    _depth = 0;
    _normal.setZero();
}

void PenaltyForceColliders::ContactPoint::initialize() {
    _colliding = false;
    _border = false;
    _processed = false;
    _tet = -1;
    _obj = -1;
    _beta.setZero();
    _depth = 0;
    _normal.setZero();
    _intersectNeighbors.clear();
    _collidingNeighbors.clear();
}

PenaltyForceColliders::IntersectionPoint::IntersectionPoint() {
    _normal.setZero();
    _t = 0;
    _depth.setZero();
}

PenaltyForceColliders::PenaltyForceColliders() = default;

PenaltyForceColliders::PenaltyForceColliders(const SoftBody &object) {
    objects.push_back(object);
    pointers.push_back(0);
    _pointList.resize(1);
    _pointList[0].resize(object.V.rows());
    _intersectingEdgeSet.resize(1);
    _exactPoints.resize(1);
    _collidingPoints.resize(1);
}

PenaltyForceColliders::PenaltyForceColliders(const std::vector<SoftBody> &colliders) {
    pointers.resize(colliders.size());
    pointers[0] = 0;
    objects.push_back(colliders[0]);
    for (size_t i = 1; i < colliders.size(); i++) {
        objects.push_back(colliders[i]);
        pointers[i] = pointers[i-1] + 3 * colliders[i-1].V.rows();
    }
    _pointList.resize(colliders.size());
    for (size_t i = 0; i < _pointList.size(); i++) {
        _pointList[i].resize(colliders[i].V.rows());
    }
    _intersectingEdgeSet.resize(colliders.size());
    _exactPoints.resize(colliders.size());
    _collidingPoints.resize(colliders.size());
}

PenaltyForceColliders::PenaltyForceColliders(
        const std::vector<SoftBody> &colliders,
        const std::vector<int> &ptrs,
        const std::vector<Boundary> &boundaries) {
    objects = colliders;
    pointers = ptrs;
    boundaryList = boundaries;
    _pointList.resize(colliders.size());
    for (size_t i = 0; i < _pointList.size(); i++) {
        _pointList[i].resize(colliders[i].V.rows());
    }
    _intersectingEdgeSet.resize(colliders.size());
    _exactPoints.resize(colliders.size());
    _collidingPoints.resize(colliders.size());
}

void PenaltyForceColliders::setStiffness(double k) {
    _stiffness = k;
}

void PenaltyForceColliders::setDamping(double d) {
    _damping = d;
}

void PenaltyForceColliders::collisionCull(const SpatialHash &hashTable, const Eigen::VectorXd& q) {
    for (int id = 0; id < objects.size(); id++) {
        const SoftBody& object = objects[id];
        Eigen::Ref<const Eigen::MatrixXi> T = object.T;
        int initialIndex = pointers[id];

        for (int t = 0; t < T.rows(); t++) {
            Eigen::RowVector4i tet = T.row(t);
            Eigen::Vector3d node0 = q.segment(initialIndex+3*tet[0], 3);
            Eigen::Vector3d node1 = q.segment(initialIndex+3*tet[1], 3);
            Eigen::Vector3d node2 = q.segment(initialIndex+3*tet[2], 3);
            Eigen::Vector3d node3 = q.segment(initialIndex+3*tet[3], 3);

            BoundingBoxD AABB(node0, node1, node2, node3);
            std::vector<size_t> keys;
            hashTable.coveredCells(keys, AABB);

            for (const size_t key : keys) {
                for (const auto &pair : hashTable[key]) {
                    if (pair.second == id) {
                        // self collision
                        if (pair.first != tet[0] && pair.first != tet[1] && pair.first != tet[2] && pair.first != tet[3]) {
                            // the point is not belong to the tet
                            Eigen::Vector3d position = q.segment(initialIndex+3*pair.first, 3);
                            if (AABB.contains(position)) {
                                Eigen::Vector4d beta = barycentricCoordinate(node0, node1, node2, node3, position);
                                if (beta[0] > 0 && beta[1] > 0 && beta[2] > 0 && beta[3] > 0) {
                                    _pointList[id][pair.first]._colliding = true;
                                    _pointList[id][pair.first]._tet = t;
                                    _pointList[id][pair.first]._obj = id;
                                    _pointList[id][pair.first]._beta = beta;
                                    _collidingPoints[id].push_back(pair.first);
                                }
                            }
                        }
                    } else {
                        // collision with other object
                        Eigen::Vector3d position = q.segment(pointers[pair.second]+3*pair.first, 3);
                        if (AABB.contains(position)) {
                            Eigen::Vector4d beta = barycentricCoordinate(node0, node1, node2, node3, position);
                            if (beta[0] > 0 && beta[1] > 0 && beta[2] > 0 && beta[3] > 0) {
                                _pointList[pair.second][pair.first]._colliding = true;
                                _pointList[pair.second][pair.first]._tet = t;
                                _pointList[pair.second][pair.first]._obj = id;
                                _pointList[pair.second][pair.first]._beta = beta;
                                _collidingPoints[pair.second].push_back(pair.first);
                            }
                        }
                    }
                }
            }
        }
    }
}

void PenaltyForceColliders::collisionForce(
        Eigen::VectorXd &f,
        const SpatialHash &hashTable,
        const Eigen::VectorXd &q,
        const Eigen::VectorXd &qdot) {
    // for multi bodies
    assemblePenaltyForce(f, qdot);
}

double PenaltyForceColliders::collisionEnergy(
        SpatialHash &hashTable,
        const Eigen::VectorXd &q,
        const Eigen::VectorXd &qdot) {
    double energy = 0;

    // for multi bodies
    // generate response energy
    energy += penaltyEnergy();

    return energy;
}

// multi body simulation
void PenaltyForceColliders::computeConstraintSet(SpatialHash& hashTable, const Eigen::VectorXd &q) {
    // cache all vertices to the hash table
    hashTable.clear();
    for (int id = 0; id < objects.size(); id++) {
        const SoftBody& obj = objects[id];
        int initialIndex = pointers[id];
        hashTable.build(q.segment(initialIndex, 3*obj.V.rows()), id);
    }
    SpatialHash vertexHash = hashTable;

    // find colliding points
    collisionCull(hashTable, q);

    // cache all intersecting edges for compute depth and normal
    constructIntersectingEdgeSet();
    hashTable.clear();
    for (int id = 0; id < objects.size(); id++) {
        const Eigen::MatrixXi& E = _intersectingEdgeSet[id];
        if (E.rows() > 0) {
            int initialIndex = pointers[id];
            const SoftBody& obj = objects[id];
            hashTable.build(E, q.segment(initialIndex, 3*obj.V.rows()), id);
        }
    }
    computeExactIntersectionPoint(hashTable, q);

    computeDepthOfBorderPoints();

    constructCollidingNeighbors();
    std::cout << "check border points: ";
    visualColliders();

    propagate(q);

    // rebuild the hash table to catch the point indices for boundaries test
    // hashTable = vertexHash;

    // for test
    // std::cout << "test hash" << std::endl;
    // hashTable.visualHashTable();
    // std::cout << "test colliding points" << std::endl;
    std::cout << "chek all colliding points: ";
    visualColliders();

}

void PenaltyForceColliders::initialize() {
    for (int id = 0; id < objects.size(); id++) {
        _collidingPoints[id].clear();
        _exactPoints[id].clear();
        _intersectingEdgeSet[id].resize(0, 2);

        for (auto &point : _pointList[id]) {
            point.initialize();
        }
    }
}

void PenaltyForceColliders::collisionStiffness(Eigen::SparseMatrixd &K, const Eigen::VectorXd &q) {
    for (int id = 0; id < objects.size(); id++) {
        int initialIndex = pointers[id];

        for (int v : _collidingPoints[id]) {
            const ContactPoint &point = _pointList[id][v];
            int collidedIndex = pointers[point._obj];
            Eigen::RowVector4i tet = objects[point._obj].T.row(point._tet);

            // compute blocks
            Eigen::Vector3d stretchSpring = point._depth * point._normal;
            // d2E / dP2
            Eigen::Matrix3d Hp = _stiffness * stretchSpring * stretchSpring.transpose() / stretchSpring.squaredNorm();
            // d2E / d0d0
            Eigen::Matrix3d H00 = point._beta[0] * Hp;
            // d2E / d1d1
            Eigen::Matrix3d H11 = point._beta[1] * Hp;
            // d2E / d2d2
            Eigen::Matrix3d H22 = point._beta[2] * Hp;
            // d2E / d3d3
            Eigen::Matrix3d H33 = point._beta[3] * Hp;

            // assemble blocks
            for (int ii = 0; ii < 3; ii++)
                for (int jj = 0; jj < 3; jj++) {
                    // diagonal elements
                    // for Hq
                    K.coeffRef(initialIndex+3*v+ii, initialIndex+3*v+jj) += Hp(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[0]+ii, collidedIndex+3*tet[0]+jj) += H00(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[1]+ii, collidedIndex+3*tet[1]+jj) += H11(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[2]+ii, collidedIndex+3*tet[2]+jj) += H22(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[3]+ii, collidedIndex+3*tet[3]+jj) += H33(ii, jj);

                    // other elements
                    K.coeffRef(initialIndex+3*v+ii, collidedIndex+3*tet[0]+jj) -= H00(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[0]+ii, initialIndex+3*v+jj) -= H00(ii, jj);

                    K.coeffRef(initialIndex+3*v+ii, collidedIndex+3*tet[1]+jj) -= H11(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[1]+ii, initialIndex+3*v+jj) -= H11(ii, jj);

                    K.coeffRef(initialIndex+3*v+ii, collidedIndex+3*tet[2]+jj) -= H22(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[2]+ii, initialIndex+3*v+jj) -= H22(ii, jj);

                    K.coeffRef(initialIndex+3*v+ii, collidedIndex+3*tet[3]+jj) -= H33(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[3]+ii, initialIndex+3*v+jj) -= H33(ii, jj);
                }
        }
    }
}

void PenaltyForceColliders::collisionDamping(Eigen::SparseMatrixd &K, const Eigen::VectorXd &qdot) {
    for (int id = 0; id < objects.size(); id++) {
        int initialIndex = pointers[id];

        for (int v: _collidingPoints[id]) {
            const ContactPoint &point = _pointList[id][v];
            int collidedIndex = pointers[point._obj];
            Eigen::RowVector4i tet = objects[point._obj].T.row(point._tet);

            // compute blocks
            Eigen::Vector3d stretchSpring = point._depth * point._normal;
            // d2E / dpdvp
            Eigen::Matrix3d Hvp = _damping * stretchSpring * stretchSpring.transpose() / stretchSpring.squaredNorm();
            // d2E / dpdv0
            Eigen::Matrix3d Hv0 = point._beta[0] * Hvp;
            // d2E / dpdv1
            Eigen::Matrix3d Hv1 = point._beta[1] * Hvp;
            // d2E / dpdv2
            Eigen::Matrix3d Hv2 = point._beta[2] * Hvp;
            // d2E / dpdv3
            Eigen::Matrix3d Hv3 = point._beta[3] * Hvp;

            // assemble blocks
            for (int ii = 0; ii < 3; ii++)
                for (int jj = 0; jj < 3; jj++) {
                    // diagonal elements
                    // for Hv
                    K.coeffRef(initialIndex+3*v+ii, initialIndex+3*v+jj) += Hvp(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[0]+ii, collidedIndex+3*tet[0]+jj) += Hv0(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[1]+ii, collidedIndex+3*tet[1]+jj) += Hv1(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[2]+ii, collidedIndex+3*tet[2]+jj) += Hv2(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[3]+ii, collidedIndex+3*tet[3]+jj) += Hv3(ii, jj);

                    // other elements
                    K.coeffRef(initialIndex+3*v+ii, collidedIndex+3*tet[0]+jj) -= Hv0(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[0]+ii, initialIndex+3*v+jj) -= Hv0(ii, jj);

                    K.coeffRef(initialIndex+3*v+ii, collidedIndex+3*tet[1]+jj) -= Hv1(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[1]+ii, initialIndex+3*v+jj) -= Hv1(ii, jj);

                    K.coeffRef(initialIndex+3*v+ii, collidedIndex+3*tet[2]+jj) -= Hv2(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[2]+ii, initialIndex+3*v+jj) -= Hv2(ii, jj);

                    K.coeffRef(initialIndex+3*v+ii, collidedIndex+3*tet[3]+jj) -= Hv3(ii, jj);
                    K.coeffRef(collidedIndex+3*tet[3]+ii, initialIndex+3*v+jj) -= Hv3(ii, jj);
                }
        }
    }
}

void PenaltyForceColliders::constructIntersectingEdgeSet() {
    for (int id = 0; id < objects.size(); id++) {
        const SoftBody& object = objects[id];
        Eigen::Ref<const Eigen::MatrixXi> E = object.E;
        std::vector<std::pair<int, int>> edges;

        for (int e = 0; e < E.rows(); e++) {
            int e0 = E(e, 0);
            int e1 = E(e, 1);

            // border point is the first end point of the intersection edge
            if (_pointList[id][e0]._colliding && !_pointList[id][e1]._colliding) {
                _pointList[id][e0]._border = true;
                _pointList[id][e0]._intersectNeighbors.push_back(edges.size());
                edges.emplace_back(e0, e1);
            }
            if (!_pointList[id][e0]._colliding && _pointList[id][e1]._colliding) {
                _pointList[id][e1]._border = true;
                _pointList[id][e1]._intersectNeighbors.push_back(edges.size());
                edges.emplace_back(e1, e0);
            }
        }

        Eigen::MatrixXi Edges;
        Edges.resize(edges.size(), 2);
        for (int i = 0; i < Edges.rows(); i++) {
            Edges(i, 0) = edges[i].first;
            Edges(i, 1) = edges[i].second;
        }
        _intersectingEdgeSet[id] = Edges;
        _exactPoints[id].resize(Edges.rows());
    }
}

void PenaltyForceColliders::constructCollidingPoints() {
    for (size_t id = 0; id < _pointList.size(); id++) {
        for (int i = 0; i < _pointList[id].size(); i++) {
            if (_pointList[id][i]._colliding) {
                _collidingPoints[id].push_back(i);
            }
        }
    }
}

void PenaltyForceColliders::constructCollidingNeighbors() {
    for (int id = 0; id < objects.size(); id++) {
        const SoftBody& object = objects[id];
        Eigen::Ref<const Eigen::MatrixXi> E = object.E;

        for (int e = 0; e < E.rows(); e++) {
            int e0 = E(e, 0);
            int e1 = E(e, 1);

            if (_pointList[id][e0]._colliding && _pointList[id][e1]._colliding) {
                _pointList[id][e0]._collidingNeighbors.push_back(e1);
                _pointList[id][e1]._collidingNeighbors.push_back(e0);
            }
        }

    }
}

void PenaltyForceColliders::computeExactIntersectionPoint(
        const SpatialHash &hashTable,
        const Eigen::VectorXd& q) {
    for (int id = 0; id < objects.size(); id++) {
        const SoftBody& object = objects[id];
        Eigen::Ref<const Eigen::MatrixXi> F = object.F;
        int initialIndex = pointers[id];

        for (int t = 0; t < F.rows(); t++) {
            Eigen::RowVector3i tri = F.row(t);
            Eigen::Vector3d node0 = q.segment(initialIndex+3*tri[0], 3);
            Eigen::Vector3d node1 = q.segment(initialIndex+3*tri[1], 3);
            Eigen::Vector3d node2 = q.segment(initialIndex+3*tri[2], 3);
            Eigen::Vector3d X20 = node2 - node0;
            Eigen::Vector3d X10 = node1 - node0;
            // should be outward
            Eigen::Vector3d localNormal = X10.cross(X20);
            localNormal.normalized();

            BoundingBoxD AABB(node0, node1, node2);
            std::vector<size_t> keys;
            hashTable.coveredCells(keys, AABB);
            for (const auto &key : keys) {
                for (const auto &pair : hashTable[key]) {
                    Eigen::RowVector2i edge = _intersectingEdgeSet[pair.second].row(pair.first);
                    // I set the first end of the intersection edge is the border point
                    int end0 = edge[0]; // border
                    int end1 = edge[1]; // non-colliding
                    if (pair.second == id) {
                        // self collision
                        if (!SoftBody::isEdgeConnectTriangle(edge, tri)) {
                            Eigen::Vector3d pos0 = q.segment(initialIndex+3*end0, 3);
                            Eigen::Vector3d pos1 = q.segment(initialIndex+3*end1, 3);
                            BoundingBoxD bboxE(pos0, pos1);
                            if (bboxE.overlaps(AABB)) {
                                // compute exact intersection point
                                Ray ray(pos1, pos0 - pos1);
                                Eigen::Vector3d beta;
                                if (ray.intersectPoint(beta, node0, node1, node2)) {
                                    IntersectionPoint& exactPoint = _exactPoints[id][pair.first];
                                    if (ray.t < 1 - exactPoint._t) {
                                        exactPoint._t = 1 - ray.t;
                                        exactPoint._depth = exactPoint._t * (pos1 - pos0);
                                        exactPoint._normal = localNormal;
                                    }
                                }
                            }
                        }
                    } else {
                        Eigen::Vector3d pos0 = q.segment(pointers[pair.second]+3*end0, 3);
                        Eigen::Vector3d pos1 = q.segment(pointers[pair.second]+3*end1, 3);
                        BoundingBoxD bboxE(pos0, pos1);
                        if (bboxE.overlaps(AABB)) {
                            // compute exact intersection point
                            Ray ray(pos1, pos0 - pos1);
                            Eigen::Vector3d beta;
                            if (ray.intersectPoint(beta, node0, node1, node2)) {
                                IntersectionPoint& exactPoint = _exactPoints[pair.second][pair.first];
                                if (ray.t < 1 - exactPoint._t) {
                                    exactPoint._t = 1 - ray.t;
                                    exactPoint._depth = exactPoint._t * (pos1 - pos0);
                                    exactPoint._normal = localNormal;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void PenaltyForceColliders::computeDepthOfBorderPoints() {
    for (int id = 0; id < objects.size(); id++) {
        for (size_t i = 0; i < _collidingPoints[id].size(); i++) {
            size_t vertex = _collidingPoints[id][i];
            if (_pointList[id][vertex]._border) {
                ContactPoint &point = _pointList[id][vertex];
                double numeratorDepth = 0;
                Eigen::Vector3d numeratorNormal;
                numeratorNormal.setZero();
                double denominator = 0;
                for (int n = 0; n < point._intersectNeighbors.size(); n++) {
                    size_t edgeNum = point._intersectNeighbors[n];
                    double distance = _exactPoints[id][edgeNum]._depth.norm();

                    double w = 1 / (distance * distance);
                    if (distance < 1e-6) {
                        w = 0;
                    }

                    // depth
                    numeratorDepth += w * _exactPoints[id][edgeNum]._depth.dot(_exactPoints[id][edgeNum]._normal);

                    // normal
                    numeratorNormal += w * _exactPoints[id][edgeNum]._normal;

                    denominator += w;
                }
                if (denominator == 0) {
                    point._depth = 0;
                    point._normal.setZero();
                } else {
                    point._depth = numeratorDepth / denominator;
                    point._normal = numeratorNormal / denominator;
                    point._normal.normalize();
                }
                point._processed = true;

                if (std::isnan(point._depth)) {
                    std::cout << "nan depth at " << vertex << std::endl;
                }
            }
        }
    }
}

void PenaltyForceColliders::propagate(const Eigen::VectorXd& q) {
    size_t numOfUnprocessedPoints;
    do {
        numOfUnprocessedPoints = 0;
        for (int id = 0; id < objects.size(); id++) {
            int initialIndex = pointers[id];
            for (size_t i = 0; i < _collidingPoints[id].size(); i++) {
                size_t vertex = _collidingPoints[id][i];
                ContactPoint &point = _pointList[id][vertex];
                if (!point._border && !point._processed) {
                    numOfUnprocessedPoints++;
                    Eigen::Vector3d pointPos = q.segment(initialIndex+3*vertex, 3);
                    double numeratorDepth = 0;
                    Eigen::Vector3d numeratorNormal;
                    double denominator = 0;
                    for (size_t n = 0; n < point._collidingNeighbors.size(); n++) {
                        size_t neighborNum = point._collidingNeighbors[n];
                        const ContactPoint &neighbor = _pointList[id][neighborNum];
                        if (neighbor._processed && neighbor._border) {
                            point._border = true;
                            Eigen::Vector3d neighborPos = q.segment(initialIndex+3*neighborNum, 3);
                            double distance = (neighborPos - pointPos).norm();
                            double mu = 1 / distance * distance;

                            // depth
                            numeratorDepth += mu * ((neighborPos - pointPos).dot(neighbor._normal) + neighbor._depth);

                            // normal
                            numeratorNormal += mu * neighbor._normal;

                            denominator += mu;
                        }
                    }
                    if (point._border) {
                        point._depth = numeratorDepth / denominator;
                        point._normal = numeratorNormal / denominator;
                        point._normal.normalized();
                        point._processed = true;
                    }

                    if (std::isnan(point._depth)) {
                        std::cout << "nan depth at " << vertex << std::endl;
                    }
                }
            }
        }
    } while (numOfUnprocessedPoints != 0);
}

void PenaltyForceColliders::assemblePenaltyForce(Eigen::VectorXd &f, const Eigen::VectorXd& qdot) {
    for (int id = 0; id < objects.size(); id++) {
        for (size_t i = 0; i < _collidingPoints[id].size(); i++) {
            size_t vertex = _collidingPoints[id][i];
            const ContactPoint &point = _pointList[id][vertex];
            Eigen::RowVector4i collidedTet = objects[point._obj].T.row(point._tet);
            int collidingIndex = pointers[id];
            int collidedIndex = pointers[point._obj];
            Eigen::Vector3d velocity0 = qdot.segment(collidedIndex+3*collidedTet[0], 3);
            Eigen::Vector3d velocity1 = qdot.segment(collidedIndex+3*collidedTet[1], 3);
            Eigen::Vector3d velocity2 = qdot.segment(collidedIndex+3*collidedTet[2], 3);
            Eigen::Vector3d velocity3 = qdot.segment(collidedIndex+3*collidedTet[3], 3);
            Eigen::Vector3d velocityOfPoint = qdot.segment(collidingIndex+3*vertex, 3);
            Eigen::Vector4d beta = point._beta;
            Eigen::Vector3d velocityOfTet = beta[0] * velocity0 + beta[1] * velocity1 + beta[2] * velocity2 + beta[3] * velocity3;
            // relative velocity from collided point(obj: A) to colliding point(obj: B) project to normal direction(from A to B)
            Eigen::Vector3d velocity = velocityOfPoint - velocityOfTet;
            double velocityInNormal = velocity.dot(point._normal);

            double magnitude = - (_stiffness * point._depth - _damping * velocityInNormal);
            // spring force of colliding point along with normal direction
            Eigen::Vector3d springForceOfPoint = - magnitude * point._normal;
            Eigen::Vector3d springForceOfTet0 = beta[0] * magnitude * point._normal;
            Eigen::Vector3d springForceOfTet1 = beta[1] * magnitude * point._normal;
            Eigen::Vector3d springForceOfTet2 = beta[2] * magnitude * point._normal;
            Eigen::Vector3d springForceOfTet3 = beta[3] * magnitude * point._normal;

            // assemble spring forces
            f.segment(collidedIndex+3*collidedTet[0], 3) += springForceOfTet0;
            f.segment(collidedIndex+3*collidedTet[1], 3) += springForceOfTet1;
            f.segment(collidedIndex+3*collidedTet[2], 3) += springForceOfTet2;
            f.segment(collidedIndex+3*collidedTet[3], 3) += springForceOfTet3;
            f.segment(collidingIndex+3*vertex, 3) += springForceOfPoint;
        }
    }
}

double PenaltyForceColliders::penaltyEnergy() {
    double energy = 0;
    for (int id = 0; id < objects.size(); id++) {
        for (size_t i = 0; i < _collidingPoints[id].size(); i++) {
            size_t vertex = _collidingPoints[id][i];
            const ContactPoint &point = _pointList[id][vertex];
            energy += _stiffness * point._depth * point._depth / 2.;
        }
    }

    return energy;
}

double PenaltyForceColliders::collideWithBoundary(
        Eigen::VectorXd& f,
        const SpatialHash &hashTable,
        const Eigen::VectorXd& q,
        const Eigen::VectorXd& qdot) {
    double energy = 0;
    for (const auto &boundary : boundaryList) {
        // DCD
        for (const auto &key : boundary.hashKeys) {
            for (const auto &pair : hashTable[key]) {
                int initialIndex = pointers[pair.second];
                Eigen::Vector3d position = q.segment(initialIndex+3*pair.first, 3);
                double signDistance = (position - boundary.point).dot(boundary.normal);
                if (signDistance < 0) {
                    Eigen::Vector3d velocity = qdot.segment(initialIndex+3*pair.first, 3);
                    double relativeV = velocity.dot(boundary.normal);
                    double distance = - signDistance;
                    double magnitude = - _stiffness * distance - _damping * relativeV;
                    Eigen::Vector3d force = magnitude * boundary.normal;
                    f.segment(initialIndex+3*pair.first, 3) += force;
                    energy += _stiffness * distance * distance / 2;
                }
            }
        }
    }

    return energy;
}

void PenaltyForceColliders::boundaryForce(Eigen::VectorXd &f, const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot) {
    f.resize(q.size());
    f.setZero();
    // DCD
    for (int id = 0; id < objects.size(); id++) {
        int initialIndex = pointers[id];
        int num = objects[id].V.rows();
        Eigen::Ref<const Eigen::VectorXd> positions = q.segment(initialIndex, 3*num);
        Eigen::Ref<const Eigen::VectorXd> velocities = qdot.segment(initialIndex, 3*num);

        /*
        // serial
        for (int n = 0; n < num; n++) {
            Eigen::Vector3d position = positions.segment(3*n, 3);
            for (const auto& boundary : boundaryList) {
                double signDistance = (position - boundary.point).dot(boundary.normal);
                if (signDistance < 0) {
                    Eigen::Vector3d velocity = velocities.segment(3*n, 3);
                    double relativeV = velocity.dot(boundary.normal);
                    double distance = -signDistance;
                    double magnitude = -(_stiffness * distance - _damping * relativeV);
                    Eigen::Vector3d force = -magnitude * boundary.normal;
                    f.segment(initialIndex + 3 * n, 3) += force;
                }
            }
        }
         */
        // concurrent
        tbb::parallel_for(tbb::blocked_range<int>(0, num), [&](const tbb::blocked_range<int> &r) {
            for (int n = r.begin(); n < r.end(); n++) {
                Eigen::Vector3d position = positions.segment(3 * n, 3);
                for (const auto &boundary: boundaryList) {
                    double signDistance = (position - boundary.point).dot(boundary.normal);
                    if (signDistance < 0) {
                        Eigen::Vector3d velocity = velocities.segment(3 * n, 3);
                        double relativeV = velocity.dot(boundary.normal);
                        double distance = -signDistance;
                        double magnitude = -(_stiffness * distance - _damping * relativeV);
                        Eigen::Vector3d force = -magnitude * boundary.normal;
                        f.segment(initialIndex + 3 * n, 3) += force;
                    }
                }
            }
        });
    }
}

double PenaltyForceColliders::boundaryEnergy(const SpatialHash& hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot) {
    double energy = 0;
    for (const auto &boundary : boundaryList) {
        // DCD
        for (const auto &key: boundary.hashKeys) {
            for (const auto &pair: hashTable[key]) {
                int initialIndex = pointers[pair.second];
                Eigen::Vector3d position = q.segment(initialIndex + 3 * pair.first, 3);
                double signDistance = (position - boundary.point).dot(boundary.normal);
                if (signDistance < 0) {
                    double distance = -signDistance;
                    energy += _stiffness * distance * distance / 2;
                }
            }
        }
    }

    return energy;
}

double PenaltyForceColliders::collisionResponse(Eigen::VectorXd &f, const SpatialHash &hashTable, const Eigen::VectorXd& q, const Eigen::VectorXd& qdot) {

    constructIntersectingEdgeSet();

    computeExactIntersectionPoint(hashTable, q);

    computeDepthOfBorderPoints();

    constructCollidingNeighbors();

    propagate(q);

    double energy = penaltyEnergy();

    assemblePenaltyForce(f, qdot);

    return energy;
}

// for test
void PenaltyForceColliders::visualColliders() const {
    for (int id = 0; id < _collidingPoints.size(); id++) {
        if (!_collidingPoints[id].empty()) {
            std::cout << "at obj " << id << " has colliding points: ";
            for (int p = 0; p < _collidingPoints[id].size(); p++) {
                int num = _collidingPoints[id][p];
                std::cout << "[number = " << num << ", " << " depth = " << _pointList[id][num]._depth << ", normal = "
                    << "("<< _pointList[id][num]._normal.x() << ", " << _pointList[id][num]._normal.y() << ", "
                    << _pointList[id][num]._normal.z() << "), " << _pointList[id][num]._normal.norm() << "], which collide with tet "
                    << _pointList[id][num]._tet << " in obj " << _pointList[id][num]._obj << "; ";
            }
            std::cout << std::endl;
        }
    }
}