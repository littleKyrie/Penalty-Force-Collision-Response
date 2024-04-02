//
// Created by Kyrie Zhang on 2023/12/7.
//

#ifndef SOFT_BODY_H
#define SOFT_BODY_H

#include <EigenTypes.h>
#include <Eigen/Dense>
#include <vector>
#include <memory>

// in the multi-body simulation,
// this class is used to conserve the state of every object in the reference space
class SoftBody {
    //data
public:
    Eigen::MatrixXd V; // vertex set
    Eigen::MatrixXi F; // triangle surface mesh
    Eigen::MatrixXi T; // tetrahedron volume mesh
    Eigen::MatrixXi E; // edge set for triangle mesh or tetrahedron mesh
    // update these items in out functions
    Eigen::VectorXd volumes;

    // constructors
public:
    SoftBody();

    SoftBody(const Eigen::MatrixXd& Vertex, const Eigen::MatrixXi& Face, const Eigen::MatrixXi& Volume, const Eigen::MatrixXi& Edge);

    SoftBody(const Eigen::MatrixXd& Vertex, const Eigen::MatrixXi& Face, const Eigen::MatrixXi& Volume);

    SoftBody(const SoftBody& other);

    // functions
public:
    // move the object
    void translate(Eigen::Ref<const Eigen::Vector3d> translation);

    void initSelfState(Eigen::VectorXd& q, Eigen::VectorXd& qdot) const;

    SoftBody& operator=(const SoftBody& other);

    // tetrahedron mesh
    void getEdgeSetForVolume();

    // triangle mesh
    void getEdgeSetForSurface();

    // determine if the edge is belonged to the triangle
    static bool isTriangleEdge(Eigen::Ref<const Eigen::RowVector2i> edge, Eigen::Ref<const Eigen::RowVector3i> tri);

    // determine if the edge is belonged to the tetrahedral
    static bool isTetrahedralEdge(Eigen::Ref<const Eigen::RowVector2i> edge, Eigen::Ref<const Eigen::RowVector4i> tet);

    // determine if the edge has connected with the triangle
    static bool isEdgeConnectTriangle(Eigen::Ref<const Eigen::RowVector2i> edge, Eigen::Ref<const Eigen::RowVector3i> tri);
};

typedef std::shared_ptr<SoftBody> SoftBodyPtr;

#endif //SOFT_BODY_H
