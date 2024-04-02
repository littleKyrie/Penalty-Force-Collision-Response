//
// Created by Kyrie Zhang on 2023/11/18.
//

#ifndef TESTSOLVER_H
#define TESTSOLVER_H

#include <EigenTypes.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <thread>
#include <oneapi/tbb.h>
#include <Gauss_Seidel_Solver.h>
#include <testTools.h>
#include <random>

void testSolver() {
    Eigen::SparseMatrixd A;
    Eigen::VectorXd b;
    A.resize(3, 3);
    A.insert(0, 0) = 8;
    A.insert(0, 1) = -3;
    A.insert(0, 2) = 2;
    A.insert(1, 0) = 4;
    A.insert(1, 1) = 11;
    A.insert(1, 2) = -1;
    A.insert(2, 0) = 6;
    A.insert(2, 1) = 3;
    A.insert(2, 2) = 12;
    b.resize(3);
    b[0] = 20;
    b[1] = 33;
    b[2] = 36;

    Eigen::SparseMatrixd AT;
    Eigen::VectorXd bT;
    AT.resize(10, 10);
    AT.insert(0, 0) = 3;
    AT.insert(0, 7) = 3;
    AT.insert(1, 1) = 15;
    AT.insert(1, 2) = 7;
    AT.insert(1, 6) = 8;
    AT.insert(2, 1) = 7;
    AT.insert(2, 2) = 28;
    AT.insert(2, 3) = 9;
    AT.insert(2, 4) = 9;
    AT.insert(2, 8) = 3;
    AT.insert(3, 2) = 9;
    AT.insert(3, 3) = 18;
    AT.insert(3, 8) = 9;
    AT.insert(4, 2) = 9;
    AT.insert(4, 4) = 20;
    AT.insert(4, 7) = 8;
    AT.insert(4, 9) = 3;
    AT.insert(5, 5) = 8;
    AT.insert(5, 9) = 8;
    AT.insert(6, 1) = 8;
    AT.insert(6, 6) = 14;
    AT.insert(6, 8) = 6;
    AT.insert(7, 0) = 3;
    AT.insert(7, 4) = 8;
    AT.insert(7, 7) = 11;
    AT.insert(8, 2) = 3;
    AT.insert(8, 3) = 9;
    AT.insert(8, 6) = 6;
    AT.insert(8, 8) = 18;
    AT.insert(9, 4) = 3;
    AT.insert(9, 5) = 8;
    AT.insert(9, 9) = 11;
    bT.resize(10);
    bT[0] = -1;
    bT[1] = 0;
    bT[2] = 5;
    bT[3] = 10;
    bT[4] = 6;
    bT[5] = 16;
    bT[6] = 7;
    bT[7] = 20;
    bT[8] = 24;
    bT[9] = 18;

    Eigen::SparseMatrixd P;
    P.resize(4, 4);
    P.insert(0, 0) = 2;
    P.insert(0, 1) = -1;
    P.insert(1, 0) = -1;
    P.insert(1, 1) = 2;
    P.insert(1, 2) = -1;
    P.insert(2, 1) = -1;
    P.insert(2, 2) = 2;
    P.insert(2, 3) = -1;
    P.insert(3, 2) = -1;
    P.insert(3, 3) = 2;
    Eigen::VectorXd r;
    r.resize(4);
    r[0] = 3;
    r[1] = -4;
    r[2] = 4;
    r[3] = -3;
    CCSd P_CCS(4, 4);
    std::vector<tbb::concurrent_vector<std::tuple<size_t, size_t, double>>> columnList1;
    columnList1.resize(4);
    for (int i = 0; i < P_CCS.rows(); i++) {
        if (i > 0) {
            columnList1[i-1].push_back(std::tuple<size_t, size_t, double>(i, i-1, -1));
            // largeP.insert(i, i-1) = -1;
        }
        if (i < P_CCS.rows()-1) {
            columnList1[i+1].push_back(std::tuple<size_t, size_t, double>(i, i+1, -1));
            // largeP.insert(i, i+1) = -1;
        }
        columnList1[i].push_back(std::tuple<size_t, size_t, double>(i, i, 2));
        //largeP.insert(i, i) = 2;
    }
    P_CCS.setFromColumnList(columnList1);
    /*
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            std::cout << P_CCS.coeff(i, j) << ", ";
        }
        std::cout << std::endl;
    }
    */

    Eigen::SparseMatrixd largeP;
    largeP.resize(699, 699);
    for (int i = 0; i < largeP.rows(); i++) {
        if (i > 0) {
            largeP.insert(i, i-1) = -1;
        }
        if (i < largeP.rows()-1) {
            largeP.insert(i, i+1) = -1;
        }
        largeP.insert(i, i) = 3;
    }
    Eigen::VectorXd largeB;
    largeB.resize(largeP.cols());
    for (size_t i = 0; i < largeB.size(); i++) {
        if (i == 0) {
            largeB[i] = 3;
        } else if (i == largeB.size()-1) {
            largeB[i] = -3;
        } else if (i % 2 == 1) {
            largeB[i] = -4;
        } else if (i % 2 == 0) {
            largeB[i] = 4;
        }
    }

    // test::visualSparseMatrix(largeP);
    // test::visualVector(largeB);

    CCSd largeP_CCS(699, 699);
    std::vector<tbb::concurrent_vector<std::tuple<size_t, size_t, double>>> columnList;
    columnList.resize(699);
    for (int i = 0; i < largeP_CCS.rows(); i++) {
        if (i > 0) {
            columnList[i-1].push_back(std::tuple<size_t, size_t, double>(i, i-1, -1));
            // largeP.insert(i, i-1) = -1;
        }
        if (i < largeP.rows()-1) {
            columnList[i+1].push_back(std::tuple<size_t, size_t, double>(i, i+1, -1));
            // largeP.insert(i, i+1) = -1;
        }
        columnList[i].push_back(std::tuple<size_t, size_t, double>(i, i, 2));
        //largeP.insert(i, i) = 2;
    }
    largeP_CCS.setFromColumnList(columnList);

    Eigen::SparseMatrixd M;
    // bool isInverseM = jacobiPreconditioning(M, tmp_H);
    // bool isInverseM = choleskyPreconditioning(M, AT);

    // a random diagonally dominant matrix for jacobi test
    int n = 699;
    Eigen::SparseMatrix<double> SpMat(n, n);
    Eigen::VectorXd rhs(n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0., 100.);

    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                SpMat.coeffRef(i, j) = dis(gen);
                sum += SpMat.coeff(i, j);
            }
        }
        SpMat.coeffRef(i, i) = dis(gen) + sum;
        rhs[i] += dis(gen);
    }

    std::chrono::high_resolution_clock::time_point start, end;
    Eigen::VectorXd dir;
    start = std::chrono::high_resolution_clock::now();
    // jacobiSolver(dir, AT, bT, false);
    jacobi(dir, SpMat, rhs, true);
    jacobiSolver(dir, SpMat, rhs);
    // gaussSeidel(dir, largeP, largeB, true);
    // conjugateGradient(dir, largeP_CCS, largeB);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> frameTime = end - start;
    // std::cout << "the time in CG of CCS version is " << frameTime.count();
    // std::cout << "the time in GS is " << frameTime.count();
    std::cout << "the time is " << frameTime.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // conjugateGradient(dir, largeP, largeB);
    end = std::chrono::high_resolution_clock::now();
    frameTime = end - start;
    // std::cout << "the time in CG of Eigen Sparse version is " << frameTime.count();
    // preconditionedConjugateGradient(dir, tmp_H, tmp_g, M, isInverseM);
    for (size_t i = 0; i < dir.size(); i++) {
        // std::cout << "x" << i << " = " << dir[i] << std::endl;
    }

    Eigen::VectorXd x;
    x.resize(b.size());
    x.setZero();
    // gaussSeidel(x, largeP, largeB, true);

    /*
    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "x" << i << " = " << x[i] << std::endl;
    }
     */
}

#endif //TESTSOLVER_H
