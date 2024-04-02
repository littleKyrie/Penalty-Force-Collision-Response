#include <Eigen/Dense>
#include <vector>
#include <colliders.h>

void find_min_vertices(
        std::vector<size_t> &indices,
        Eigen::Ref<const Eigen::MatrixXd> V,
        double tol = 1e-3);

void findMultiMinVertices(
        std::vector<std::vector<size_t>>& indices,
        const std::vector<SoftBody>& objects,
        const std::vector<int>& fixedObjs,
        double tol = 1e-3);