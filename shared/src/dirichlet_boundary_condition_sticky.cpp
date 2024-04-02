//
// Created by Kyrie Zhang on 2024/3/7.
//

#include <dirichlet_boundary_condition_stiky.h>
#include <oneapi/tbb.h>

void readFixedVertices(std::vector<size_t>& bcNodes, const std::string file) {
    std::cout << "reading dirichlet boundary condition nodes from " << file << std::endl;

    std::ifstream bcNodeIn;
    bcNodeIn.open(file);
    if (!bcNodeIn.good()) {
        std::cout << "can not open the file " << file << std::endl;
        return;
    }

    // deal with the first line to get the total number of bc nodes
    int numOfBCs = 0;
    std::string line;
    std::getline(bcNodeIn, line);
    std::stringstream ss(line);
    char c;
    ss >> c;
    ss >> numOfBCs;

    while (1) {
        std::getline(bcNodeIn, line);
        // check the end of the input stream
        if (bcNodeIn.eof()) {
            break;
        }
        // skip the empty line
        if (line.empty()) {
            continue;
        }
        // skip the comment
        if (line.at(0) == '#') {
            continue;
        }

        std::stringstream sstream(line);
        int bcNode;
        sstream >> bcNode;
        bcNodes.push_back(bcNode);
    }

    bcNodeIn.close();
}

void findTopLeftVertices(std::vector<size_t>& bcNodes, const Eigen::MatrixXd& V, double tol) {
    // first find top nodes
    std::vector<size_t> topNodes;
    double maxVertex = V(0, 1);

    for (int i = 0; i < V.rows(); i++) {
        if (maxVertex < V(i, 1)) {
            maxVertex = V(i, 1);
        }
    }

    for (int i = 0; i < V.rows(); i++) {
        if (std::abs(maxVertex - V(i, 1)) <= tol) {
            topNodes.push_back(i);
        }
    }

    // second proceed to find left nodes
    double leftmostMaxVertex = V(topNodes[0], 0);
    for (size_t topNode : topNodes) {
        if (leftmostMaxVertex > V(topNode, 0)) {
            leftmostMaxVertex = V(topNode, 0);
        }
    }

    for (size_t topNode : topNodes) {
        if (std::abs(leftmostMaxVertex - V(topNode, 0)) <= tol) {
            bcNodes.push_back(topNode);
        }
    }
}

void findTopBackVertices(std::vector<size_t>& bcNodes, const Eigen::MatrixXd& V, double tol) {
    // first find top nodes
    std::vector<size_t> topNodes;
    double maxVertex = V(0, 1);

    for (int i = 0; i < V.rows(); i++) {
        if (maxVertex < V(i, 1)) {
            maxVertex = V(i, 1);
        }
    }

    for (int i = 0; i < V.rows(); i++) {
        if (std::abs(maxVertex - V(i, 1)) <= tol) {
            topNodes.push_back(i);
        }
    }

    // second proceed to find back nodes
    double leftmostMaxVertex = V(topNodes[0], 2);
    for (size_t topNode : topNodes) {
        if (leftmostMaxVertex > V(topNode, 2)) {
            leftmostMaxVertex = V(topNode, 2);
        }
    }

    for (size_t topNode : topNodes) {
        if (std::abs(leftmostMaxVertex - V(topNode, 2)) <= tol) {
            bcNodes.push_back(topNode);
        }
    }
}

void findBackVertices(std::vector<size_t>& bcNodes, const Eigen::MatrixXd& V, double tol) {
    double backNode = V(0, 2);

    for (int i = 0; i < V.rows(); i++) {
        if (backNode > V(i, 2)) {
            backNode  = V(i, 2);
        }
    }

    for (int i = 0; i < V.rows(); i++) {
        if (std::abs(backNode - V(i, 2)) < tol) {
            bcNodes.push_back(i);
        }
    }
}

void findForwardVertices(std::vector<size_t>& bcNodes, const Eigen::MatrixXd& V, double tol) {
    double forwardNode = V(0, 2);

    for (int i = 0; i < V.rows(); i++) {
        if (forwardNode < V(i, 2)) {
            forwardNode = V(i, 2);
        }
    }

    for (int i = 0; i < V.rows(); i++) {
        if (std::abs(forwardNode - V(i, 2)) < tol) {
            bcNodes.push_back(i);
        }
    }
}

void constraintSize(Eigen::SparseMatrixd &mat, const Eigen::SparseMatrixd &tmpMat, const std::vector<int> &map) {
    Eigen::SparseMatrixd sparseMat(tmpMat);
    int n = map.size();

    tbb::parallel_for(tbb::blocked_range<int>(0, n*n), [&](const tbb::blocked_range<int> &r) {
        for (int i = r.begin(); i < r.end(); i++) {
            int vi = i / n;
            int vj = i % n;

            int row = map[vi];
            int col = map[vj];

            if (row != -1 && col != -1) {
                for (size_t jj = 0; jj < 3; jj++) {
                    for (Eigen::SparseMatrixd::InnerIterator it(sparseMat, 3 * col + jj); it; ++it) {
                        if (it.row() == 3 * row) {
                            it.valueRef() = mat.coeff(vi, vj + jj);
                            ++it;
                            it.valueRef() = mat.coeff(vi + 1, vj + jj);
                            ++it;
                            it.valueRef() = mat.coeff(vi + 2, vj + jj);
                            break;
                        }
                    }
                }
            }
        }
    });

    // write back
    mat = sparseMat;
}