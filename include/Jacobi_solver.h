//
// Created by Kyrie Zhang on 2023/10/20.
// some solver for linear equation system

#ifndef JACOBI_SOLVER_H
#define JACOBI_SOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

#include <iostream>
#include <thread>
#include <vector>

#include <oneapi/tbb.h>
#include <oneapi/tbb/enumerable_thread_specific.h>

template <typename T>
bool jacobi(
        Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
        const Eigen::SparseMatrix<T>& A,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        bool isConcurrent = true,
        size_t maxSteps = 200,
        double tol = 1e-4,
        double alpha = 0.3) {
    // initialize
    size_t size = b.size();
    x.resize(size); // final result
    x.setZero();

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double, std::milli> Time = end - start;

    start = std::chrono::high_resolution_clock::now();
    // spectral radius of iterative matrix to determine if the solver will converge
    // compute iterative matrix B
    Eigen::SparseMatrix<T> D_inv(b.size(), b.size());
    D_inv.setIdentity();
    if (!isConcurrent) {
        for (size_t i = 0; i < b.size(); ++i) {
            D_inv.coeffRef(i, i) = 1 / A.coeff(i, i);
        }
    } else {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, b.size()), [&](const tbb::blocked_range<size_t>& r){
            for (size_t i = r.begin(); i < r.end(); i++) {
                D_inv.coeffRef(i, i) = 1 / A.coeff(i, i);
            }
        });
    }
    Eigen::SparseMatrix<T> E(b.size(), b.size());
    E.setIdentity();

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B = E - alpha * A * D_inv;
    Eigen::EigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> eigenSolver(B);
    T spectralRadius = std::numeric_limits<T>::min();
    for (size_t i = 0; i < eigenSolver.eigenvalues().size(); i++) {
        // std::cout << "the eigen value at " << i << "= " << eigenvalues[i] << std::endl;
        if (abs(eigenSolver.eigenvalues()[i]) > spectralRadius) {
            spectralRadius = abs(eigenSolver.eigenvalues()[i]);
        }
    }
    std::cout << spectralRadius << std::endl;
    if (spectralRadius >= 1) {
        std::cout << "the Jacobi can not converge" << std::endl;
        std::cout << spectralRadius << std::endl;
        return false;
    }
    end = std::chrono::high_resolution_clock::now();
    Time = end - start;
    // std::cout << "time of computing the spectral radius = " << Time.count() << std::endl;

    // jacobi solver
    // initialize variables
    Eigen::Matrix<T, Eigen::Dynamic, 1> y(size); // intermediate variable
    y.setZero(); // initial guess
    Eigen::Matrix<T, Eigen::Dynamic, 1> r; // residual error
    r.resize(size);
    r.setZero();

    // for chebyshev acceleration
    double w = 0;
    // reserve y(k-1)
    Eigen::Matrix<T, Eigen::Dynamic, 1> last_y = y;

    start = std::chrono::high_resolution_clock::now();
    size_t iterator = 0;
    while (iterator < maxSteps) {
        // reserve y(k)
        Eigen::Matrix<T, Eigen::Dynamic, 1> old_y = y;

        // update y(k) to obtain y(k+1)
        if (!isConcurrent) {
            for (size_t i = 0; i < size; i++) {
                double sum = 0;
                for (size_t j = 0; j < size; j++) {
                    sum += A.coeff(i, j) * old_y[j];
                }
                r[i] = b[i] - sum;

                // update the i-th components with relaxation
                y[i] = old_y[i] + alpha * r[i] / A.coeff(i, i);
            }
        } else {
            tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const tbb::blocked_range<size_t>& k){
                for (size_t i = k.begin(); i < k.end(); i++) {
                    double sum = 0;
                    for (size_t j = 0; j < size; j++) {
                        sum += A.coeff(i, j) * old_y[j];
                    }
                    r[i] = b[i] - sum;

                    // update the i-th components with relaxation
                    y[i] = old_y[i] + alpha * r[i] / A.coeff(i, i);
                }
            });
        }

        // determine converge with r(k)
        if (r.norm() < tol) {
            break;
        }

        // chebyshev acceleration
        if (iterator == 0) {
            w = 1;
        } else if (iterator == 1) {
            w = 2 / (2 - spectralRadius * spectralRadius);
        } else {
            w = 4 / (4 - spectralRadius * spectralRadius * w);
        }

        y = w * y + (1 - w) * last_y;
        last_y = old_y;

        iterator++;

        // std::cout << "at iteration " << iterator << std::endl;
        // std::cout << "the residual is " << r.norm() << std::endl;
    }
    end = std::chrono::high_resolution_clock::now();
    Time = end - start;
    // std::cout << "time of solving the system = " << Time.count();

    x = y;

    std::cout << "at iteration " << iterator << std::endl;
    std::cout << "the residual is " << r.norm() << std::endl;

    return true;
}

// may be ignored
template <typename T>
bool jacobiSolver (
        Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
        const Eigen::SparseMatrix<T>& A,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        bool isConcurrent = true,
        T alpha = 0.8,
        size_t maxSteps = 200,
        T tol = 1e-3) {
    // initialize
    size_t size = b.size();
    x.resize(size); // final result
    x.setZero();

    // spectral radius of iterative matrix to determine if the solver will converge
    // compute iterative matrix B
    Eigen::SparseMatrix<T> D_inv(b.size(), b.size());
    D_inv.setIdentity();
    if (!isConcurrent) {
        for (size_t i = 0; i < b.size(); ++i) {
            D_inv.coeffRef(i, i) = 1 / A.coeff(i, i);
        }
    } else {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, b.size()), [&](const tbb::blocked_range<size_t>& r){
            for (size_t i = r.begin(); i < r.end(); i++) {
                D_inv.coeffRef(i, i) = 1 / A.coeff(i, i);
            }
        });
    }
    Eigen::SparseMatrix<T> E(b.size(), b.size());
    E.setIdentity();

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B = E - alpha * A * D_inv;
    Eigen::EigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> eigenSolver(B);
    T spectralRadius = std::numeric_limits<T>::min();
    for (size_t i = 0; i < eigenSolver.eigenvalues().size(); i++) {
        // std::cout << "the eigen value at " << i << "= " << eigenvalues[i] << std::endl;
        if (abs(eigenSolver.eigenvalues()[i]) > spectralRadius) {
            spectralRadius = abs(eigenSolver.eigenvalues()[i]);
        }
    }
    std::cout << spectralRadius << std::endl;
    if (spectralRadius >= 1) {
        std::cout << "the Jacobi can not converge" << std::endl;
        std::cout << spectralRadius << std::endl;
        return false;
    }

    // jacobi solver
    // initialize variables
    Eigen::Matrix<T, Eigen::Dynamic, 1> y(size); // intermediate variable
    y.setZero(); // initial guess
    Eigen::Matrix<T, Eigen::Dynamic, 1> r; // residual error
    r.resize(size);
    r.setZero();

    // for chebyshev acceleration
    double w = 0;
    // reserve y(k-1)
    Eigen::Matrix<T, Eigen::Dynamic, 1> last_y = y;

    size_t iterator = 0;
    while (iterator < maxSteps) {
        // reserve y(k)
        Eigen::Matrix<T, Eigen::Dynamic, 1> old_y = y;

        // update y(k) to obtain y(k+1)
        if (!isConcurrent) {
            for (size_t i = 0; i < size; i++) {
                double sum = 0;
                for (size_t j = 0; j < size; j++) {
                    sum += A.coeff(i, j) * old_y[j];
                }
                r[i] = b[i] - sum;

                // update the i-th components with relaxation
                y[i] = old_y[i] + alpha * r[i] / A.coeff(i, i);
            }
        } else {
            tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const tbb::blocked_range<size_t>& k){
                for (size_t i = k.begin(); i < k.end(); i++) {
                    double sum = 0;
                    for (size_t j = 0; j < size; j++) {
                        sum += A.coeff(i, j) * old_y[j];
                    }
                    r[i] = b[i] - sum;

                    // update the i-th components with relaxation
                    y[i] = old_y[i] + alpha * r[i] / A.coeff(i, i);
                }
            });
        }

        // determine converge with r(k)
        if (r.norm() < tol) {
            break;
        }

        // chebyshev acceleration
        if (iterator == 0) {
            w = 1;
        } else if (iterator == 1) {
            w = 2 / (2 - spectralRadius * spectralRadius);
        } else {
            w = 4 / (4 - spectralRadius * spectralRadius * w);
        }

        // y = w * y + (1 - w) * last_y;
        // last_y = old_y;

        iterator++;

        // std::cout << "at iteration " << iterator << std::endl;
        // std::cout << "the residual is " << r.norm() << std::endl;
    }

    x = y;

    std::cout << "at iteration " << iterator << std::endl;
    std::cout << "the residual is " << r.norm() << std::endl;

    return true;
}

#endif //JACOBI_SOLVER_H
