//
// Created by Kyrie Zhang on 2023/11/10.
//

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

#include <EigenTypes.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <iostream>
#include <thread>

#include <oneapi/tbb.h>

template <typename dataType>
bool jacobiPreconditioning(
        Eigen::SparseMatrix<dataType>& M_inv,
        const Eigen::SparseMatrix<dataType>& A) {
    M_inv.resize(A.rows(), A.cols());
    M_inv.setZero();

    for (size_t i = 0; i < A.rows(); i++) {
        M_inv.insert(i, i) = 1 / A.coeff(i, i);
    }

    return true; // return the inverse M
}

template <typename dataType>
bool SGSPreconditioning(
        Eigen::SparseMatrix<dataType>& M,
        const Eigen::SparseMatrix<dataType>& A,
        dataType w = 1) {


    return false;
}

template <typename dataType>
bool choleskyPreconditioning(
        Eigen::SparseMatrix<dataType>& M,
        const Eigen::SparseMatrix<dataType>& A) {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<dataType>> cholesky;
    cholesky.analyzePattern(A);
    cholesky.factorize(A);

    Eigen::SparseMatrix<dataType> L;
    // Eigen::SparseMatrix<dataType> U;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> d;
    if (cholesky.info() == Eigen::Success) {
        L = cholesky.matrixL();
        // U = cholesky.matrixU();
        d = cholesky.vectorD();
    } else {
        std::cout << "fail to factorize A by LDLT" << std::endl;
        exit(0);
    }

    tbb::parallel_for(tbb::blocked_range<size_t>(0, d.size()),
            [&](const tbb::blocked_range<size_t>& r){
                for (size_t i = r.begin(); i < r.end(); i++) {
                    d[i] = std::sqrt(d[i]);
                }
    });
    // Eigen::Matrix<dataType, Eigen::Dynamic, Eigen::Dynamic> D = d.asDiagonal();
    Eigen::SparseMatrix<dataType> DT;
    DT.resize(d.size(), d.size());

    /*
    tbb::parallel_for(tbb::blocked_range<size_t>(0, d.size()),
            [&](const tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i < r.end(); i++) {
            DT.coeffRef(i, i) = d[i];
        }
    });
     */

    for (size_t i = 0; i < d.size(); i++) {
        DT.insert(i, i) = d[i];
    }

    M = L * DT;

    /*
    std::cout << "-----matrix M-----" << std::endl;
    for (size_t i = 0; i < d.size(); i++) {
        for (size_t j = 0; j < d.size(); j++) {
            std::cout << M.coeff(i, j) << " ";
        }
        std::cout << std::endl;
    }
     */

    return false; // return LD^1/2, not its inverse
}

#endif //PRECONDITIONER_H
