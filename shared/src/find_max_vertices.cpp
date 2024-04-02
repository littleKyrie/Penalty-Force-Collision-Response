//
// Created by Kyrie Zhang on 2023/12/16.
//

#include <find_max_vertices.h>

void find_max_vertices(
        std::vector<size_t> &indices,
        Eigen::Ref<const Eigen::MatrixXd> V,
        double tol) {

    double max_vertex = V(0,1);

    for(unsigned int vi=0; vi<V.rows(); ++vi) {
        max_vertex = (V(vi,1) > max_vertex ? V(vi,1) : max_vertex);
    }

    for(unsigned int vi=0; vi<V.rows(); ++vi) {

        if(std::abs(V(vi,1)-max_vertex) <= tol) {
            indices.push_back(vi);
        }
    }
}

void findMultiMaxVertices(
        std::vector<std::vector<size_t>>& indices,
        const std::vector<SoftBody>& objects,
        const std::vector<int>& fixedObjs,
        double tol) {
    indices.resize(fixedObjs.size());
    for (int n = 0; n < fixedObjs.size(); n++) {
        int id = fixedObjs[n];
        const SoftBody& obj = objects[id];
        Eigen::Ref<const Eigen::MatrixXd> V = obj.V;

        double maxVertex = V(0, 1);

        for(unsigned int vi=0; vi<V.rows(); ++vi) {
            maxVertex = (V(vi,1) > maxVertex ? V(vi,1) : maxVertex);
        }

        for(unsigned int vi=0; vi<V.rows(); ++vi) {
            if(std::abs(V(vi,1)-maxVertex) <= tol) {
                indices[n].push_back(vi);
            }
        }
    }
}