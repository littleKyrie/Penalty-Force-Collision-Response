//
// Created by Kyrie Zhang on 2023/11/15.
//

#ifndef TESTCCS_H
#define TESTCCS_H

#include <EigenTypes.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <SparseStorage.h>

#include <iostream>
#include <vector>
#include <oneapi/tbb.h>
#include <thread>

/*
 * use a template sparse matrix as same with the eigen document
 *       0   3   0   0   0
 *       22  0   0   0   17
 *       7   5   0   1   0
 *       0   0   0   0   0
 *       0   0   14  0   8
 */

void testEigenSparse () {
    Eigen::SparseMatrix<double> E(5, 5);
    E.insert(1, 0) = 0;
    E.insert(2, 0) = 0;
    E.insert(0, 1) = 3;
    E.insert(2, 1) = 0;
    E.insert(4, 2) = 0;
    E.insert(2, 3) = 0;
    E.insert(1, 4) = 0;
    E.insert(4, 4) = 0;
    E.makeCompressed();

    std::chrono::high_resolution_clock::time_point s, e;
    s = std::chrono::high_resolution_clock::now();
    int* outer = E.outerIndexPtr();
    int* inner = E.innerIndexPtr();
    double* value = E.valuePtr();
    e = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> localTime = e - s;
    std::cout << "get time of create ptr: " << localTime.count() << std::endl;

    s = std::chrono::high_resolution_clock::now();
    Eigen::SparseMatrixd::InnerIterator it(E, 0);
    e = std::chrono::high_resolution_clock::now();
    localTime = e - s;
    std::cout << "get time of InnerIterator: " << localTime.count() << std::endl;

    std::cout << "the matrix E is " << std::endl << Eigen::MatrixXd(E) << std::endl;
    std::cout << "it has " << E.nonZeros() << " none zero elements " << std::endl;

    std::vector<double> blocks;
    blocks.push_back(22);
    blocks.push_back(7);
    blocks.push_back(3);
    blocks.push_back(5);
    blocks.push_back(14);
    blocks.push_back(1);
    blocks.push_back(17);
    blocks.push_back(8);

    // E.resize(5, 5);
    // E.setZero();
    Eigen::Ref<Eigen::SparseMatrixd> Z = E;
    // Z.coeffRef(2,0) = 7;
    // std::cout << "the matrix E is " << std::endl << Eigen::MatrixXd(E) << std::endl;
    // std::cout << "it has " << E.nonZeros() << " none zero elements " << std::endl;

    int index = 0;
    for (Eigen::SparseMatrixd::InnerIterator it(E, 1); it; ++it) {
        for (size_t i = 0; i < 2; i++) {
            if (it.row() == 2) {
                it.valueRef() = it.value() + 7;
                break;
            }
            ++index;
        }
    }
    std::cout << "the matrix E is " << std::endl << Eigen::MatrixXd(E) << std::endl;
    std::cout << "it has " << E.nonZeros() << " none zero elements " << std::endl;
    std::cout << index;
}

void testColumnList() {
    std::vector<tbb::concurrent_vector<std::tuple<size_t, size_t, double>>> columnList;
    CCSd A(5, 5);
    columnList.resize(5);
    std::tuple<size_t, size_t, double> triplet0(1, 0, 22);
    std::tuple<size_t, size_t, double> triplet1(2, 0, 7);
    std::tuple<size_t, size_t, double> triplet2(0, 1, 3);
    std::tuple<size_t, size_t, double> triplet3(2, 1, 5);
    std::tuple<size_t, size_t, double> triplet4(4, 2, 14);
    std::tuple<size_t, size_t, double> triplet5(2, 3, 1);
    std::tuple<size_t, size_t, double> triplet6(1, 4, 17);
    std::tuple<size_t, size_t, double> triplet7(4, 4, 8);
    columnList[0].push_back(triplet1);
    columnList[0].push_back(triplet0);
    columnList[1].push_back(triplet3);
    columnList[1].push_back(triplet2);
    columnList[2].push_back(triplet4);
    columnList[3].push_back(triplet5);
    columnList[4].push_back(triplet6);
    columnList[4].push_back(triplet7);
    /*
    std::cout << "unordered list" << std::endl;
    for (size_t i = 0; i < columnList.size(); i++) {
        for (const auto &ele : columnList[i]) {
            std::cout << "(" << std::get<0>(ele) << ", " << std::get<1>(ele) << ", " << std::get<2>(ele) << "), ";
        }
        std::cout << std::endl;
    }
     */

    /*
    for (size_t i = 0; i < columnList.size(); i++) {
        std::sort(columnList[i].begin(), columnList[i].end(),
                  [](const std::tuple<size_t, size_t, double>& ele1,
                          const std::tuple<size_t, size_t, double>& ele2){
           return std::get<0>(ele1) < std::get<0>(ele2);
        });
    }
    std::cout << "ordered list" << std::endl;
    for (size_t i = 0; i < columnList.size(); i++) {
        for (const auto &ele : columnList[i]) {
            std::cout << "(" << std::get<0>(ele) << ", " << std::get<1>(ele) << ", " << std::get<2>(ele) << "), ";
        }
        std::cout << std::endl;
    }
     */

    A.setFromColumnList(columnList);
    std::cout << "show the matrix A" << std::endl;
    for (size_t i = 0; i < A.rows(); i++) {
        for (size_t j = 0; j < A.cols(); j++) {
            std::cout << A.coeff(i, j) << ", ";
        }
        std::cout << std::endl;
    }

    /*
    A.insertValue(0, 0, 3);
    std::cout << A.coeff(0, 0) << std::endl;
    A.insertValue(0, 2, 0);
    std::cout << A.coeff(0, 2) << std::endl;
    A.insertValue(1, 0, 5);
     */
    std::vector<size_t> z(10);
    std::vector<size_t> y(10);
    for (size_t i = 0; i < 10; i++) {
        z[i] = i+1;
        y[i] = 0;
    }
    struct Scan {
        size_t sum;
        std::vector<size_t>& y;
        const std::vector<size_t>& z;

        // Scan(std::vector<size_t>& y_, const std::vector<size_t>& z_)
    };
    size_t res = 0;
    size_t result = tbb::parallel_scan(tbb::blocked_range<size_t>(0, 10), res,[&](const tbb::blocked_range<size_t>& r, size_t tmp, bool is_final_scan) -> size_t {
            size_t sum = tmp;
            for (size_t i = r.begin(); i != r.end(); ++i) {
                sum += z[i];
                if (is_final_scan) {
                    y[i] = sum;
                }
            }
            return sum;
            },
            [](size_t left, size_t right) {
                return left + right;
            }
    );
    std::cout << result << std::endl;
    for (size_t i = 0; i < 10; i++) {
        std::cout << y[i] << ", ";
    }
    std::cout << std::endl;

    CCSd B(5, 5);
    std::vector<tbb::concurrent_vector<double>> columnList2;
    columnList2.resize(5*5);
    columnList2[1].push_back(11);
    columnList2[1].push_back(11);
    columnList2[2].push_back(3);
    columnList2[2].push_back(4);
    columnList2[5].push_back(3);
    columnList2[7].push_back(2);
    columnList2[7].push_back(1);
    columnList2[7].push_back(2);
    columnList2[14].push_back(14);
    columnList2[17].push_back(1);
    columnList2[21].push_back(17);
    columnList2[24].push_back(8);
    B.setFromColumnClusters(columnList2);
    std::cout << "show the matrix B" << std::endl;
    for (size_t i = 0; i < B.rows(); i++) {
        for (size_t j = 0; j < B.cols(); j++) {
            std::cout << B.coeff(i, j) << ", ";
        }
        std::cout << std::endl;
    }
}

#endif //TESTCCS_H
