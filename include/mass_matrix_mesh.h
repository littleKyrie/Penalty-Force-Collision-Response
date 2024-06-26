#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <colliders.h>

//Input:
//  qdot - generalized velocity for the FEM system
//  T - the mx4 vertex indices for tet mesh
//  density - density of material
//  v0 - the undeformed tetrahedra volumes
//Output:
//  M - Sparse mass matrix for the whole mesh.
void mass_matrix_mesh(
        Eigen::SparseMatrixd &M,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        Eigen::Ref<const Eigen::MatrixXi> T,
        double density, Eigen::Ref<const Eigen::VectorXd> v0);

void massMatrixMultiBodies(
        Eigen::SparseMatrixd &M,
        Eigen::Ref<const Eigen::VectorXd> qdot,
        const std::vector<SoftBody>& objects,
        const std::vector<int>& pointers, double density);