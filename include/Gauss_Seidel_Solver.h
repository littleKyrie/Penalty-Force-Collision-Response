//
// Created by Kyrie Zhang on 2023/10/23.
//

#ifndef GAUSS_SEIDEL_SOLVER_H
#define GAUSS_SEIDEL_SOLVER_H

#include <EigenTypes.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <thread>
#include <oneapi/tbb.h>
#include <iostream>
#include <SparseStorage.h>
#include <graph.h>

template <typename dataType>
void gaussSeidel(Eigen::Matrix<dataType, Eigen::Dynamic, 1>& result,
                 const Eigen::SparseMatrix<dataType>& A,
                 const Eigen::Matrix<dataType, Eigen::Dynamic, 1>& b,
                 bool isConcurrent = true,
                 dataType relaxation = 1,
                 size_t maxSteps = 200,
                 dataType tol = 1e-3) {
    size_t size = b.size();
    result.resize(size);
    result.setZero();

    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<double, std::milli> Time = end - start;

    start = std::chrono::high_resolution_clock::now();
    // for graph color
    std::vector<Vertex> vertexSet;
    std::vector<Edge> edgeSet;
    std::vector<tbb::concurrent_vector<int>> concurrentSet;
    if (isConcurrent) {
        vertexSet.resize(size);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size),
                [&](const tbb::blocked_range<size_t>& r){
            for (size_t i = r.begin(); i < r.end(); i++) {
                Vertex v(i);
                vertexSet[i] = v;
            }
        });

        tbb::concurrent_vector<Edge> tmpEdgeSet;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size),
                [&](const tbb::blocked_range<size_t>& r){
            for (size_t j = r.begin(); j < r.end(); j++) {
                for (size_t i = 0; i < j; i++) {
                    if (A.coeff(i, j) != 0) {
                        Vertex v1(i);
                        Vertex v2(j);
                        Edge e(v1, v2);
                        tmpEdgeSet.push_back(e);
                    }
                }
            }
        });
        edgeSet.resize(tmpEdgeSet.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, tmpEdgeSet.size()),
                [&](const tbb::blocked_range<size_t>& r){
            for (size_t i = r.begin(); i < r.end(); i++) {
                edgeSet[i] = tmpEdgeSet[i];
            }
        });

        Graph graphColor;
        graphColor.buildGraph(vertexSet, edgeSet);
        size_t numColor = graphColor.broadFirstSearch();

        concurrentSet.resize(numColor);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size),
                [&](const tbb::blocked_range<size_t>& r){
            for (size_t i = r.begin(); i < r.end(); i++) {
                int currentColor = graphColor.visitVertexColor(i);
                int number = graphColor.visitVertexNumber(i);
                concurrentSet[currentColor].push_back(number);
            }
        });
        end = std::chrono::high_resolution_clock::now();
        Time = end - start;
        // std::cout << "the time of create grpah and color the variables = " << Time.count() << std::endl;

        // for test color
        /*
        for (size_t i = 0; i < concurrentSet.size(); i++) {
            std::cout << "at color " << i << ":" << std::endl;
            for (const auto &element : concurrentSet[i]) {
                std::cout << element << ", ";
            }
            std::cout << std::endl;
        }
         */
    }

    start = std::chrono::high_resolution_clock::now();
    // solve
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> x;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> oldX;
    oldX.resize(size);
    oldX.setZero();
    x.resize(size);
    x.setZero();

    // initial guess
    /*
    for (size_t i = 0; i < size; i++) {
        oldX[i] = 2;
        x[i] = 2;
    }
     */

    dataType residual = 0;
    Eigen::Matrix<dataType, Eigen::Dynamic, 1> resVec;
    resVec.resize(size);
    resVec.setZero();

    size_t iterator = 0;
    while (iterator < maxSteps) {
        if (!isConcurrent) {
            for (int i = 0; i < size; i++) {
                dataType sum = 0;

                for (int j = 0; j < size; j++) {
                    if (j < i) {
                        sum += A.coeff(i, j) * x[j];
                    } else if (j > i) {
                        sum += A.coeff(i, j) * oldX[j];
                    }
                }

                x[i] = (b[i] - sum) / A.coeff(i, i);

                // for SOR
                x[i] = relaxation * x[i] + (1 - relaxation) * oldX[i];
            }
        } else {
            for (const auto &set : concurrentSet) {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, set.size()),
                        [&](const tbb::blocked_range<size_t>& r){
                    for (size_t k = r.begin(); k < r.end(); k++) {
                        int i = set[k];
                        dataType sum = 0;

                        for (int j = 0; j < size; j++) {
                            if (i != j) {
                                sum += A.coeff(i, j) * x[j];
                            }
                        }

                        x[i] = (b[i] - sum) / A.coeff(i, i);

                        // for SOR
                        x[i] = relaxation * x[i] + (1 - relaxation) * oldX[i];
                    }
                });
            }
        }

        resVec = b - A * x;
        residual = resVec.norm();
        if (residual < tol)
        {
            break;
        }

        oldX = x;

        iterator++;

        // std::cout << "the GS solver iteration at " << iterator << std::endl;
        // std::cout << "and the residual of iteration is " << residual << std::endl;
    }
    end = std::chrono::high_resolution_clock::now();
    Time = end - start;
    // std::cout << "the time of soling system = " << Time.count() << std::endl;

    std::cout << "iterator is " << iterator << " of Gauss-Seidel ";
    std::cout << "and the residual is " << residual << std::endl;

    result = x;
}

#endif //GAUSS_SEIDEL_SOLVER_H
