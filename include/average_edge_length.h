//
// Created by Kyrie Zhang on 2023/12/19.
//

#ifndef AVERAGE_EDGE_LENGTH_H
#define AVERAGE_EDGE_LENGTH_H

#include <colliders.h>
#include <EigenTypes.h>
#include <vector>

double aveEdgeLength(const std::vector<SoftBody>& objects);

void getEdgeSetForVolume(Eigen::MatrixXi& E, const Eigen::MatrixXi& T);

double longestEdge(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E);

#endif //AVERAGE_EDGE_LENGTH_H
