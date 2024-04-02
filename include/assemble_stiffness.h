#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <d2V_linear_tetrahedron_dq2.h>
#include <SparseStorage.h>
#include <colliders.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  K - the sparse, global stiffness matrix
void assemble_stiffness(
        Eigen::SparseMatrixd &K,
        Eigen::Ref<const Eigen::MatrixXi> Map,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::MatrixXi> T,
        Eigen::Ref<const Eigen::VectorXd> v0,
        double C, double D);

void assembleMultiStiffness(Eigen::SparseMatrixd& K, Eigen::Ref<const Eigen::VectorXd> q, const CollidersPtr colliders, double C, double D);

// parallel assemble the elastic Hessian
void assembleHessianFast(
        Eigen::SparseMatrixd &K,
        Eigen::Ref<const Eigen::VectorXd> q,
        Eigen::Ref<const Eigen::MatrixXd> V,
        Eigen::Ref<const Eigen::MatrixXi> T,
        Eigen::Ref<const Eigen::VectorXd> v0,
        std::vector<std::map<int, Eigen::Matrix3d>>& blockList,
        double C, double D);