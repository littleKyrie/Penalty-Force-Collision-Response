//
// Created by Kyrie Zhang on 2023/12/21.
//

#ifndef TESTTOOLS_H
#define TESTTOOLS_H

#include <EigenTypes.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

namespace test {
    void visualSparseMatrix(const Eigen::SparseMatrixd& mat);

    void visualDenseMatrix(const Eigen::MatrixXd& mat);

    void visualDenseMatrix(const Eigen::MatrixXi& mat);

    void visualVector(const Eigen::VectorXd& vec);

    void visualVectorEle(const Eigen::VectorXd& vec, int ele);

    bool isMatSymmetric(const Eigen::SparseMatrixd& mat);

    bool isMatSPD(const Eigen::SparseMatrixd& mat);

    // check if there are some tetrahedrons with negative volumes
    // there are still some bugs
    void visualTetVolume(const Eigen::VectorXd &q, const Eigen::MatrixXi &T);
}

#endif //TESTTOOLS_H
