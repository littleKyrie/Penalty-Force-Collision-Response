//
// Created by Kyrie Zhang on 2023/11/30.
//

#include <init_elasticity_stiffness.h>

void initElasticityStiffness(
        Eigen::SparseMatrixd& K,
        Eigen::MatrixXi& Map,
        Eigen::Ref<const Eigen::MatrixXi> T,
        int size) {
    K.resize(size, size);

    typedef Eigen::Triplet<double> Tuple;
    std::vector<Tuple> tripletList;

    size_t numTet = T.rows();
    for (int n = 0; n < numTet; n++) {
        Eigen::Vector4i element = T.row(n);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int rowGlobal = 3 * element[i];
                int colGlobal = 3 * element[j];

                for (int ii = 0; ii < 3; ii++)
                    for (int jj = 0; jj < 3; jj++) {
                        tripletList.emplace_back(rowGlobal+ii, colGlobal+jj, 0);
                    }
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    K.makeCompressed();

    Map.resize(size, size);
    for (size_t i = 0; i < Map.rows(); i++) {
        for (size_t j = 0; j < Map.cols(); j++) {
            Map(i, j) = -1;
        }
    }
    int* outer = K.outerIndexPtr();
    int* inner = K.innerIndexPtr();
    for (int j = 0; j < Map.cols(); j++) {
        int colIndex = j;
        int start = outer[colIndex];
        int end = outer[colIndex+1];

        for (int i = start; i < end; i++) {
            int rowIndex = inner[i];
            // reserve the relative index in the column
            Map(rowIndex, j) = i - start;
        }
    }
}

void initElasticityStiffness(
        Eigen::SparseMatrixd &K,
        std::vector<int> &map,
        const std::vector<size_t> &bcNodes,
        Eigen::Ref<const Eigen::MatrixXi> T,
        int size) {
    // set the new size of K and record the number of unconstrained nodes
    int m = size - 3 * bcNodes.size();
    int n = m / 3;
    K.resize(m, m);
    // map is used to preserve the new positions of the unconstrained blocks in K
    map.resize(n);
    for (int i = 0; i < n; i++) {
        map[i] = 0;
    }

    // prepare for filling the blocks with 0 matrix
    typedef Eigen::Triplet<double> Tuple;
    std::vector<Tuple> tripletList;

    for (size_t node : bcNodes) {
        map[node] = -1;
    }
    int k = -1; // preserve the last unconstrained node's position
    for (int i = 0; i < n; i++) {
        if (map[i] == 0) {
            map[i] = ++k;
        } else if (map[i] == -1) {
            continue;
        }
    }

    // fill blocks
    for (int t = 0; t < T.rows(); t++) {
        Eigen::Vector4i element = T.row(t);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int vi = 3 * element[i];
                int vj = 3 * element[j];

                int row = map[vi];
                int col = map[vj];

                if (row != -1 && col != -1) {
                    for (int ii = 0; ii < 3; ii++)
                        for (int jj = 0; jj < 3; jj++) {
                            tripletList.emplace_back(row+ii, col+jj, 0);
                        }
                }
            }
        }
    }

    K.setFromTriplets(tripletList.begin(), tripletList.end());
    K.makeCompressed();
}

void initMultiBodiesStiffness(Eigen::SparseMatrixd& K, const std::vector<SoftBody>& objects, const std::vector<int>& pointers, int size) {
    typedef Eigen::Triplet<double> Tuple;
    std::vector<Tuple> tripletList;

    K.resize(size, size);

    for (int n = 0; n < objects.size(); n++) {
        const SoftBody& obj = objects[n];
        int initialIndex = pointers[n];
        Eigen::Ref<const Eigen::MatrixXi> T = obj.T;

        for (int t = 0; t < T.rows(); t++) {
            Eigen::RowVector4i element = T.row(t);

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++) {
                    int rowGlobal = initialIndex + 3 * element[i];
                    int colGlobal = initialIndex + 3 * element[j];

                    for (int ii = 0; ii < 3; ii++)
                        for (int jj = 0; jj < 3; jj++) {
                            tripletList.emplace_back(rowGlobal + ii, colGlobal + jj, 0);
                        }
                }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    K.makeCompressed();
}