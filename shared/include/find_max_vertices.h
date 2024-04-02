//
// Created by Kyrie Zhang on 2023/12/16.
//

#ifndef FIND_MAX_VERTICES_H
#define FIND_MAX_VERTICES_H

#include <vector>
#include <EigenTypes.h>
#include <colliders.h>

void find_max_vertices(
        std::vector<size_t> &indices,
        Eigen::Ref<const Eigen::MatrixXd> V,
        double tol = 1e-3);

void findMaxVertices(std::vector<size_t> &indices,
                     Eigen::Ref<const Eigen::MatrixXd> V,
                     int obj, double tol = 1e-3);

void findMultiMaxVertices(
        std::vector<std::vector<size_t>>& indices,
        const std::vector<SoftBody>& objects,
        const std::vector<int>& fixedObjs,
        double tol = 1e-3);

#endif //FIND_MAX_VERTICES_H
