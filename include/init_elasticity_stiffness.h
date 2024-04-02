//
// Created by Kyrie Zhang on 2023/11/30.
//

#ifndef INIT_ELASTICITY_STIFFNESS_H
#define INIT_ELASTICITY_STIFFNESS_H

#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <colliders.h>

// for concurrently assembling complete stiffness matrix
void initElasticityStiffness(
        Eigen::SparseMatrixd& K,
        Eigen::MatrixXi& Map,
        Eigen::Ref<const Eigen::MatrixXi> T,
        int size);

// for preserving the culled matrix model constrained by the Dirichlet boundary conditions
void initElasticityStiffness(
        Eigen::SparseMatrixd& K,
        std::vector<int> &map,
        const std::vector<size_t> &bcNodes,
        Eigen::Ref<const Eigen::MatrixXi> T,
        int size);

// for concurrently assembling complete matrix of all objects
void initMultiBodiesStiffness(
        Eigen::SparseMatrixd& K,
        const std::vector<SoftBody>& objects,
        const std::vector<int>& pointers,
        int size);

#endif //INIT_ELASTICITY_STIFFNESS_H
