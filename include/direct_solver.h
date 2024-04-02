//
// Created by Kyrie Zhang on 2023/11/10.
//

#ifndef DIRECT_SOLVER_H
#define DIRECT_SOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

#include <iostream>
#include <thread>

#include <oneapi/tbb.h>

// for computing the inverse matrix used in preconditioning CG
template <typename dataType>
bool gaussElimination(
        Eigen::Matrix<dataType, Eigen::Dynamic, 1>& x,
        const Eigen::SparseMatrix<dataType>& L,
        const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& b) {
    const Eigen::SparseMatrix<dataType>& U = L.transpose();
    assert(L.rows() == L.cols()
        && U.rows() == U.cols()
        && L.rows() == U.rows()
        && L.rows() == b.size());

    size_t size = b.size();
    x.resize(size);
    x.setZero();

    // forward elimination
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> y;
    y.resize(size);
    y.setZero();
    for (size_t i = 0; i < size; i++) {
        if (L.coeff(i, i) == 0) {
            std::cout << "the forward matrix is singular " << std::endl;
            return false;
        } else {
            dataType sum = 0;
            for (size_t j = 0; j < i; j++) {
                sum += L.coeff(i, j) * y[j];
            }
            y[i] = (b[i] - sum) / L.coeff(i, i);
        }
    }

    // backward substitution
    int n = int(size) - 1;
    for (int i = n; i >= 0; i--) {
        if (U.coeff(i, i) == 0) {
            std::cout << "the backward matrix is singular " << std::endl;
            return false;
        } else {
            dataType sum = 0;
            for (int j = n; j > i; j--) {
                sum += U.coeff(i, j) * x[j];
            }
            x[i] = (y[i] - sum) / U.coeff(i, i);
        }
    }

    return true; // success solved
}

// this is the new version of the up function and the real name of this algorithm is back and forward solver
template <typename T>
bool backForwardSolver(
        Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
        const Eigen::SparseMatrix<T>& L,
        const Eigen::SparseMatrix<T>& U,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) {
    int n = b.size();
    x.resize(n);
    x.setZero();
    Eigen::Matrix<T, Eigen::Dynamic, 1> y(n);
    y.setZero();

    // solve Ly = b
    for (int i = 0; i < n; i++) {
        T sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L.coeff(i, j) * y[j];
        }
        y[i] = (b[i] - sum) / L.coeff(i, i);
    }

    // solve Ux = y
    for (int i = n - 1; i >= 0; i--) {
        T sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U.coeff(i, j) * x[j];
        }
        x[i] = (y[i] - sum) / U.coeff(i, i);
    }

    return true;
}

#endif //DIRECT_SOLVER_H
