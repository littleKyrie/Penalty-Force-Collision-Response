//
// Created by Kyrie Zhang on 2023/11/9.
//

#ifndef CONJUGATE_GRADIENT_SOLVER_H
#define CONJUGATE_GRADIENT_SOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

#include <iostream>
#include <vector>
#include <thread>

#include <direct_solver.h>
#include <SparseStorage.h>

template<typename dataType>
bool conjugateGradient(
        Eigen::Matrix<dataType, Eigen::Dynamic, 1>& result,
        const Eigen::SparseMatrix<dataType>& A,
        const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& b,
        dataType tol = 1e-20) {
    size_t maxSteps = b.size();
    // size_t maxSteps = 10;

    result.resize(b.size());
    result.setZero();
    // Eigen::Matrix<dataType, Eigen::Dynamic, 1> x0 = result;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> x = result;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> res;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> dir;
    res.resize(b.size());
    res = b - A * x;
    dir.resize(b.size());
    dir = res;
    size_t iterator = 0;
    dataType alpha = 0;
    dataType beta = 0;
    dataType resNorm = res.norm();
    dataType res2Norm = res.squaredNorm();

    while (iterator < maxSteps) {
        if (resNorm < tol) {
            break;
        }

        Eigen::Matrix<dataType, Eigen::Dynamic, 1> conjugatedDir = A * dir;
        alpha = res2Norm / (dir.transpose() * conjugatedDir);
        x += alpha * dir;

        if (std::isnan(x.norm())) {
            std::cout << "the CG solver has numerical issue" << std::endl;
            return false;
        }

        iterator++;
        if (iterator % 50 == 0) {
            res = b - A * x;
        } else {
            res -= alpha * conjugatedDir;
        }

        dataType newRes2Norm = res.squaredNorm();

        beta = newRes2Norm / res2Norm;
        dir = res + beta * dir;

        res2Norm = newRes2Norm;
        resNorm = res.norm();
    }

    std::cout << "the number of iteration is " << iterator << std::endl;
    if (iterator == maxSteps && resNorm > tol) {
        std::cout << "the norm of residual in CG is " << resNorm << std::endl;
        std::cout << "the CG solver may not converge" << std::endl;
        // return false;
    }

    result = x;

    return true;
}

template <typename dataType>
bool preconditionedConjugateGradient(
        Eigen::Matrix<dataType, Eigen::Dynamic, 1>& result,
        const Eigen::SparseMatrix<dataType>& A,
        const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& b,
        const Eigen::SparseMatrix<dataType>& M,
        bool isInverseM,
        dataType tol = 1e-3) {
    size_t maxSteps = b.size();

    result.resize(b.size());
    result.setZero();
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> x = result;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> res;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> dir;
    res.resize(b.size());
    res = b - A * x;
    dir.resize(b.size());
    dir.setZero();
    if (isInverseM) {
        dir = M * res;
    } else {
        bool solved = gaussElimination(dir, M, res);
        if (!solved) {
            std::cout << "some issue in gauss elimination of computing initial director" << std::endl;
            return false;
        }
    }
    size_t iterator = 0;
    dataType alpha = 0;
    dataType beta = 0;
    dataType deltaNew = res.transpose() * dir;
    dataType delta0 = deltaNew;
    // dataType resNorm = res.norm();
    // dataType res2Norm = res.squaredNorm();

    // use relative residual error to determine if the system has converged
    dataType rhsNorm = b.norm();
    dataType relativeRes = 0;
    relativeRes = res.norm() / rhsNorm;

    while (iterator < maxSteps) {
        /*
        if (deltaNew <= tol * tol * delta0) {
            break;
        }
         */
        if (relativeRes <= tol) {
            break;
        }

        Eigen::Matrix<dataType, Eigen::Dynamic, 1> q = A * dir;
        alpha = deltaNew / (dir.transpose() * q);
        x += alpha * dir;

        if (std::isnan(x.norm())) {
            std::cout << "the CG solver has numerical issue" << std::endl;
            return false;
        }

        iterator++;
        if (iterator % 50 == 0) {
            res = b - A * x;
        } else {
            res -= alpha * q;
        }

        Eigen::Matrix<dataType, Eigen::Dynamic, 1> s;
        s.resize(b.size());
        s.setZero();
        if (isInverseM) {
            s = M * res;
        } else {
            bool solved = gaussElimination(s, M, res);
            if (!solved) {
                std::cout << "some issue in gauss elimination in iteration" << std::endl;
                return false;
            }
        }

        dataType deltaOld = deltaNew;
        deltaNew = res.transpose() * s;

        beta = deltaNew / deltaOld;
        dir = s + beta * dir;

        relativeRes = res.norm() / rhsNorm;
    }

    // std::cout << "the number of iteration is " << iterator << std::endl;
    if (iterator == maxSteps && relativeRes > tol) {
        // std::cout << "the norm of residual in CG is " << std::sqrt(deltaNew) << std::endl;
        std::cout << "the norm of residual in CG is " << relativeRes << std::endl;
        std::cout << "the CG solver may not converge" << std::endl;
        // return false;
    }

    result = x;

    return true;
}

// compute M_inv implicitly in this function by Eigen
template <typename T>
bool PCG(Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
         const Eigen::SparseMatrix<T>& A,
         const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
         const std::string& preconditioner = "Jacobi",
         T tol = 1e-3) {
    // set size of the problem
    size_t n = b.size();

    // initialize result
    x.resize(n);
    x.setZero();

    // prepare universal data
    size_t i = 0;
    Eigen::Matrix<T, Eigen::Dynamic, 1> x0(n);
    x0.setZero();
    Eigen::Matrix<T, Eigen::Dynamic, 1> r = b - A * x0;
    Eigen::Matrix<T, Eigen::Dynamic, 1> d(n);
    d.setZero();
    T rhsNorm = b.norm();
    T relativeResidual = r.norm() / rhsNorm;

    if (preconditioner == "Jacobi" || preconditioner == "diagonal") {
        typedef Eigen::Triplet<T> triplet;
        // diagonal preconditioner
        Eigen::SparseMatrix<T> D(n, n);
        std::vector<triplet> diagCoefficient;
        for (int k = 0; k < n; k++) {
            diagCoefficient.emplace_back(k, k, 1 / A.coeff(k, k));
        }
        D.setFromTriplets(diagCoefficient.begin(), diagCoefficient.end());

        // set initial search direction
        d = D * r;
        T ro = r.dot(d);
        while (i < n && relativeResidual > tol) {
            Eigen::Matrix<T, Eigen::Dynamic, 1> q = A * d;
            // compute step size along current search direction
            T alpha = ro / d.dot(q);
            // update x_k
            x0 += alpha * d;

            // update search direction
            if (i % 50 == 0) {
                r = b - A * x0;
            } else {
                r = r - alpha * q;
            }
            relativeResidual = r.norm() / rhsNorm;
            Eigen::Matrix<T, Eigen::Dynamic, 1> s = D * r;
            T ro_old = ro;
            ro = r.dot(s);
            T beta = ro / ro_old;
            d = s + beta * d;
            i++;
        }

    } else if (preconditioner == "incompleteCholesky" || preconditioner == "Cholesky") {
        // incomplete Cholesky factorization
        Eigen::IncompleteCholesky<T> ichol;
        ichol.compute(A);

        // set initial search direction
        d = ichol.solve(r);
        T ro = r.dot(d);
        while (i < n && relativeResidual > tol) {
            Eigen::Matrix<T, Eigen::Dynamic, 1> q = A * d;
            // compute step size along current search direction
            T alpha = ro / d.dot(q);
            // update x_k
            x0 += alpha * d;

            // update search direction
            if (i % 50 == 0) {
                r = b - A * x0;
            } else {
                r = r - alpha * q;
            }
            relativeResidual = r.norm() / rhsNorm;
            Eigen::Matrix<T, Eigen::Dynamic, 1> s = ichol.solve(r);
            T ro_old = ro;
            ro = r.dot(s);
            T beta = ro / ro_old;
            d = s + beta * d;
            i++;
        }
    } else if (preconditioner == "SSOR") {
        // w = 1 in default
        // set preconditioner
        Eigen::SparseMatrix<T> E(n, n);
        Eigen::SparseMatrix<T> F(n, n);
        Eigen::SparseMatrix<T> L(n, n);
        Eigen::SparseMatrix<T> U(n, n);
        Eigen::SparseMatrix<T> D(n, n);
        Eigen::SparseMatrix<T> D_inv(n, n);
        Eigen::SparseMatrix<T> I(n, n);
        I.setIdentity();

        typedef Eigen::Triplet<T> triplet;
        std::vector<triplet> lowerCoefficient;
        std::vector<triplet> upperCoefficient;
        std::vector<triplet> diagCoefficient;
        std::vector<triplet> diagInvCoefficient;
        for (int k = 0; k < n; k++) {
            for (int m = 0; m < k; m++) {
                if (A.coeff(k, m) != 0) {
                    lowerCoefficient.emplace_back(k, m, A.coeff(k, m));
                }
            }
            for (int j = k + 1; j < n; j++) {
                if (A.coeff(k, j) != 0) {
                    upperCoefficient.emplace_back(k, j, A.coeff(k, j));
                }
            }
            diagCoefficient.emplace_back(k, k, A.coeff(k, k));
            diagInvCoefficient.emplace_back(k, k, 1 / A.coeff(k, k));
        }
        E.setFromTriplets(lowerCoefficient.begin(), lowerCoefficient.end());
        F.setFromTriplets(upperCoefficient.begin(), upperCoefficient.end());
        D.setFromTriplets(diagCoefficient.begin(), diagCoefficient.end());
        D_inv.setFromTriplets(diagInvCoefficient.begin(), diagInvCoefficient.end());

        Eigen::SparseMatrix<T> M(n, n);
        M = (D + E) * D_inv * (D + F);
        // L = I + E * D_inv;
        // U = D + F;

        Eigen::SparseLU<Eigen::SparseMatrix<T>> lu;
        lu.compute(M);
        // set initial search direction
        // backForwardSolver(d, L, U, r);
        d = lu.solve(r);
        T ro = r.dot(d);
        while (i < n && relativeResidual > tol) {
            Eigen::Matrix<T, Eigen::Dynamic, 1> q = A * d;
            // compute step size along current search direction
            T alpha = ro / d.dot(q);
            // update x_k
            x0 += alpha * d;

            // update search direction
            if (i % 50 == 0) {
                r = b - A * x0;
            } else {
                r = r - alpha * q;
            }
            relativeResidual = r.norm() / rhsNorm;
            Eigen::Matrix<T, Eigen::Dynamic, 1> s;
            // backForwardSolver(s, L, U, r);
            s = lu.solve(r);
            T ro_old = ro;
            ro = r.dot(s);
            T beta = ro / ro_old;
            d = s + beta * d;
            i++;
        }
    }

    // std::cout << "iterate = " << i << "with relative residual = " << relativeResidual << std::endl;

    x = x0;
    return true;
}

template <typename T>
bool conjugateGradient(
        Eigen::Matrix<T, Eigen::Dynamic, 1>& result,
        const CompressedColumnStorage<T> & A,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        T tol = 1e-20) {
    size_t maxSteps = b.size();
    // size_t maxSteps = 10;

    result.resize(b.size());
    result.setZero();
    // Eigen::Matrix<dataType, Eigen::Dynamic, 1> x0 = result;
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = result;
    Eigen::Matrix<T, Eigen::Dynamic, 1> res;
    Eigen::Matrix<T, Eigen::Dynamic, 1> dir;
    res.resize(b.size());
    // res = b - A * x;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, A.rows()), [&](const tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i < r.end(); i++) {
            T sum = 0;
            for (size_t j = 0; j < A.cols(); j++) {
                sum += A.coeff(i, j) * x[j];
            }
            res[i] = b[i] - sum;
        }
    });
    // std::cout << "res0 = ";
    // for (size_t i = 0; i < res.size(); i++) {
        // std::cout << res[i] << ", ";
    // }
    // std::cout << std::endl;
    /*
    for (size_t i = 0; i < A.rows(); i++) {
        T sum = 0;
        for (size_t j = 0; j < A.cols(); j++) {
            sum += A.coeff(i, j) * x[j];
        }
        res[i] = b[i] - sum;
    }
     */

    dir.resize(b.size());
    dir = res;

    size_t iterator = 0;
    T alpha = 0;
    T beta = 0;
    T resNorm = res.norm();
    T res2Norm = res.squaredNorm();

    while (iterator < maxSteps) {
        // std::cout << "at " << iterator << std::endl;

        if (resNorm < tol) {
            break;
        }

        // Eigen::Matrix<T, Eigen::Dynamic, 1> conjugatedDir = A * dir;
        Eigen::Matrix<T, Eigen::Dynamic, 1> conjugatedDir;
        conjugatedDir.resize(b.size());
        conjugatedDir.setZero();
        tbb::parallel_for(tbb::blocked_range<size_t>(0, A.rows()), [&](const tbb::blocked_range<size_t>& r){
            for (size_t i = r.begin(); i < r.end(); i++) {
                for (size_t j = 0; j < A.cols(); j++) {
                    conjugatedDir[i] += A.coeff(i, j) * dir[j];
                }
            }
        });
        /*
        for (size_t i = 0; i < A.rows(); i++)
            for (size_t j = 0; j < A.cols(); j++) {
                conjugatedDir[i] += A.coeff(i, j) * dir[j];
            }
        */

        alpha = res2Norm / (dir.transpose() * conjugatedDir);

        x += alpha * dir;

        if (std::isnan(x.norm())) {
            std::cout << "the CG solver has numerical issue" << std::endl;
            return false;
        }

        iterator++;
        if (iterator % 50 == 0) {
            // res = b - A * x;
            /*
            for (size_t i = 0; i < A.rows(); i++) {
                T sum = 0;
                for (size_t j = 0; j < A.cols(); j++) {
                    sum += A.coeff(i, j) * x[j];
                }
                res[i] = b[i] - sum;
            }
             */
            tbb::parallel_for(tbb::blocked_range<size_t>(0, A.rows()), [&](const tbb::blocked_range<size_t>& r){
                for (size_t i = r.begin(); i < r.end(); i++) {
                    T sum = 0;
                    for (size_t j = 0; j < A.cols(); j++) {
                        sum += A.coeff(i, j) * x[j];
                    }
                    res[i] = b[i] - sum;
                }
            });
        } else {
            res -= alpha * conjugatedDir;
        }

        T newRes2Norm = res.squaredNorm();

        beta = newRes2Norm / res2Norm;
        dir = res + beta * dir;

        res2Norm = newRes2Norm;
        resNorm = res.norm();
    }

    std::cout << "the number of iteration is " << iterator << std::endl;
    if (iterator == maxSteps && resNorm > tol) {
        std::cout << "the norm of residual in CG is " << resNorm << std::endl;
        std::cout << "the CG solver may not converge" << std::endl;
        // return false;
    }

    result = x;

    return true;
}

#endif //CONJUGATE_GRADIENT_SOLVER_H
