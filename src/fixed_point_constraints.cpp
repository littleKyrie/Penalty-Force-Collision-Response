#include <fixed_point_constraints.h>
#include <algorithm>
#include <iostream>

void fixed_point_constraints(
        Eigen::SparseMatrixd &P,
        unsigned int q_size,
        const std::vector<size_t>& indices) {
    // there are 2 ways to construct P:
    // the first one is a 3nx3n diagonal matrix which has zero as some elements to present the position of fixed vertex
    // the other one is a 3(n-m)x3n matrix which has the same structure with selection matrix E to select non-fixed vertex
    // we choose the second one as the P in this function
    size_t numberOfFixedVertex = indices.size();
    size_t numberOfVertex = q_size / 3;
    auto n = static_cast<Eigen::Index>(numberOfVertex);
    auto m = static_cast<Eigen::Index>(numberOfFixedVertex);

    if (numberOfFixedVertex == 0) {
        return;
    }

    P.resize(3*(n-m), q_size);
    P.setZero();

    size_t k = 0; // from 0 to (n-m)
    for (size_t i = 0; i < numberOfVertex; i++) {
        bool non_fixed = true;

        for (size_t j = 0; j < numberOfFixedVertex; j++) {
            if (i == indices[j]) {
                non_fixed = false;
                break;
            }
        }

        if (non_fixed) {
            auto rowIndex = static_cast<Eigen::Index>(k);
            auto vertexIndex = static_cast<Eigen::Index>(i);
            P.insert(rowIndex*3+0, vertexIndex*3+0) = 1;
            P.insert(rowIndex*3+1, vertexIndex*3+1) = 1;
            P.insert(rowIndex*3+2, vertexIndex*3+2) = 1;

            k++;
        }
    }

    if (k != (n-m)) {
        std::cout << 3*(n - m);
        std::cout << "k = " << k << "is not equal to the correct number" << std::endl;
        exit(0);
    }

    std::cout << "finish construct constraints" << std::endl;
}

void fixedMultiPointsConstraints(
        Eigen::SparseMatrixd &P,
        const std::vector<std::vector<size_t>>& indices,
        const std::vector<int>& fixedObjs,
        const std::vector<SoftBody>& objects,
        const std::vector<int>& pointers) {
    Eigen::SparseMatrixd E = P;

    size_t numOfZeroRows = 0;
    for (int n = 0; n < fixedObjs.size(); n++) {
        int id = fixedObjs[n];
        const SoftBody& obj = objects[id];
        Eigen::Ref<const Eigen::MatrixXd> V = obj.V;
        int initialIndex = pointers[id];

        for (int v = 0; v < V.rows(); v++) {
            for (size_t index : indices[n]) {
                if (v == index) {
                    int rowGlobal = initialIndex + 3 * v;
                    int colGlobal = initialIndex + 3 * v;

                    E.coeffRef(rowGlobal, colGlobal) = 0;
                    E.coeffRef(rowGlobal+1, colGlobal+1) = 0;
                    E.coeffRef(rowGlobal+2, colGlobal+2) = 0;
                    numOfZeroRows += 3;
                }
            }
        }
    }

    P.resize(E.rows()-numOfZeroRows, E.cols());
    P.setZero();
    int row = 0;
    for (int i = 0; i < E.rows(); i++) {
        if (E.coeff(i, i) != 0) {
            P.coeffRef(row, i) = 1;
            row++;
        }
    }
}