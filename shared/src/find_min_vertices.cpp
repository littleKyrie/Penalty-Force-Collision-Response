#include <find_min_vertices.h>
#include <iostream>

void find_min_vertices(
        std::vector<size_t> &indices,
        Eigen::Ref<const Eigen::MatrixXd> V,
        double tol) {

    double min_vertex = V(0,1); 

    for(unsigned int vi=0; vi<V.rows(); ++vi) {
        min_vertex = (V(vi,1) < min_vertex ? V(vi,1) : min_vertex);
    }

    for(unsigned int vi=0; vi<V.rows(); ++vi) {

        if(std::abs(V(vi,1)-min_vertex) <= tol) {
            indices.push_back(vi);
        }
    }

    // concurrent
}

void findMultiMinVertices(std::vector<std::vector<size_t>>& indices, const std::vector<SoftBody>& objects, const std::vector<int>& fixedObjs, double tol) {
    indices.resize(fixedObjs.size());
    for (int n = 0; n < fixedObjs.size(); n++) {
        int id = fixedObjs[n];
        const SoftBody& obj = objects[id];
        Eigen::Ref<const Eigen::MatrixXd> V = obj.V;

        double minVertex = V(0, 1);

        for(unsigned int vi=0; vi<V.rows(); ++vi) {
            minVertex = (V(vi,1) < minVertex ? V(vi,1) : minVertex);
        }

        for(unsigned int vi=0; vi<V.rows(); ++vi) {
            if(std::abs(V(vi,1)-minVertex) <= tol) {
                indices[n].push_back(vi);
            }
        }
    }
}