//
// Created by Kyrie Zhang on 2023/12/21.
//

#include <testTools.h>

namespace test {
    void visualSparseMatrix(const Eigen::SparseMatrixd& mat) {
        for (int i = 0; i < mat.rows(); i++) {
            for (int j = 0; j < mat.cols(); j++) {
                std::cout << mat.coeff(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    void visualDenseMatrix(const Eigen::MatrixXd& mat) {
        for (int i = 0; i < mat.rows(); i++) {
            for (int j = 0; j < mat.cols(); j++) {
                std::cout << mat.coeff(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    void visualDenseMatrix(const Eigen::MatrixXi& mat) {
        for (int i = 0; i < mat.rows(); i++) {
            for (int j = 0; j < mat.cols(); j++) {
                std::cout << mat.coeff(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    void visualVector(const Eigen::VectorXd& vec) {
        for (int i = 0; i < vec.size(); i++) {
            std::cout << vec[i] << " ";
        }
        std::cout << std::endl;
    }

    void visualVectorEle(const Eigen::VectorXd& vec, int ele) {
        for (int i = 0; i < vec.size() / 3; i++) {
            std::cout << vec[3*i+ele] << " ";
        }
        std::cout << std::endl;
    }

    bool isMatSymmetric(const Eigen::SparseMatrixd& mat) {
        if (mat.isApprox(mat.transpose())) {
            return true;
        }

        return false;
    }

    bool isMatSPD(const Eigen::SparseMatrixd& mat) {
        if (isMatSymmetric(mat)) {
            Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrixd> solver(mat);
            if (solver.eigenvalues().minCoeff() >= 0) {
                return true;
            } else {
                return false;
            }
        }

        return false;
    }

    void visualTetVolume(const Eigen::VectorXd &q, const Eigen::MatrixXi &T) {
        int num = 0;
        for (int t = 0; t < T.rows(); t++) {
            Eigen::Vector3d x0 = q.segment(T(t, 0), 3);
            Eigen::Vector3d x1 = q.segment(T(t, 1), 3);
            Eigen::Vector3d x2 = q.segment(T(t, 2), 3);
            Eigen::Vector3d x3 = q.segment(T(t, 3), 3);

            double v = (x3 - x0).dot((x1 - x0).cross(x2 - x0)) / 6;

            if (v < 0) {
                std::cout << "at tet " << t << " : " << v << " ";
                num++;
            }
        }
        if (num > 0) {
            std::cout << std::endl;
        }
    }
}
