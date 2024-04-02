//
// Created by Kyrie Zhang on 2024/3/7.
//

#ifndef DIRICHLET_BOUNDARY_CONDITION_STIKY_H
#define DIRICHLET_BOUNDARY_CONDITION_STIKY_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <Eigen/Dense>
#include <EigenTypes.h>

void readFixedVertices(std::vector<size_t>& bcNodes, const std::string file);

// first find the top nodes and second find the left nodes on the top
void findTopLeftVertices(std::vector<size_t>& bcNodes, const Eigen::MatrixXd& V, double tol = 1e-3);

// first find the top nodes and second find the back nodes on the top
void findTopBackVertices(std::vector<size_t>& bcNodes, const Eigen::MatrixXd& V, double tol = 1e-3);

// just find the back nodes
void findBackVertices(std::vector<size_t>& bcNodes, const Eigen::MatrixXd& V, double tol = 1e-3);

// just find the forward nodes
void findForwardVertices(std::vector<size_t>& bcNodes, const Eigen::MatrixXd& V, double tol = 1e-3);

// parallel regulate the size of Matrix based on the constraints
void constraintSize(Eigen::SparseMatrixd &mat, const Eigen::SparseMatrixd &tmpMat, const std::vector<int>& map);

#endif //DIRICHLET_BOUNDARY_CONDITION_STIKY_H
