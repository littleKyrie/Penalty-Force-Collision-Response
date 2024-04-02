//
// Created by Kyrie Zhang on 2023/12/19.
//

#include <average_edge_length.h>

double aveEdgeLength(const std::vector<SoftBody>& objects) {
    double length = 0;
    int numOfEdges = 0;

    for (int id = 0; id < objects.size(); id++) {
        const SoftBody& obj = objects[id];
        Eigen::Ref<const Eigen::MatrixXi> E = obj.E;
        Eigen::Ref<const Eigen::MatrixXd> V = obj.V;

        for (int e = 0; e < E.rows(); e++) {
            Eigen::RowVector3d p0 = V.row(E(e, 0));
            Eigen::RowVector3d p1 = V.row(E(e, 1));

            double localLength = (p0 - p1).norm();
            length += localLength;
            numOfEdges++;
        }
    }

    length = length / numOfEdges;
    return length;
}

void getEdgeSetForVolume(Eigen::MatrixXi& E, const Eigen::MatrixXi& T) {
    std::vector<std::pair<int, int>> tmp_edge;
    for (int t = 0; t < T.rows(); t++) {
        Eigen::RowVector4i tet = T.row(t);
        int p0 = tet[0];
        int p1 = tet[1];
        int p2 = tet[2];
        int p3 = tet[3];

        if (p0 < p1) {
            tmp_edge.emplace_back(p0, p1);
        } else {
            tmp_edge.emplace_back(p1, p0);
        }
        if (p0 < p2) {
            tmp_edge.emplace_back(p0, p2);
        } else {
            tmp_edge.emplace_back(p2, p0);
        }
        if (p0 < p3) {
            tmp_edge.emplace_back(p0, p3);
        } else {
            tmp_edge.emplace_back(p3, p0);
        }
        if (p1 < p2) {
            tmp_edge.emplace_back(p1, p2);
        } else {
            tmp_edge.emplace_back(p2, p1);
        }
        if (p1 < p3) {
            tmp_edge.emplace_back(p1, p3);
        } else {
            tmp_edge.emplace_back(p3, p1);
        }
        if (p2 < p3) {
            tmp_edge.emplace_back(p2, p3);
        } else {
            tmp_edge.emplace_back(p3, p2);
        }
    }
    std::sort(tmp_edge.begin(), tmp_edge.end(), [&](const std::pair<int, int>& e0, const std::pair<int, int>& e1) -> bool {
        if (e0.first < e1.first) {
            return true;
        } else if (e0.first == e1.first && e0.second < e1.second) {
            return true;
        }
        return false;
    });

    std::vector<std::pair<int, int>> edgeSet;
    edgeSet.push_back(tmp_edge[0]);
    for (size_t i = 1; i < tmp_edge.size(); i++) {
        if (tmp_edge[i-1].first != tmp_edge[i].first || tmp_edge[i-1].second != tmp_edge[i].second) {
            edgeSet.push_back(tmp_edge[i]);
        }
    }

    E.resize(edgeSet.size(), 2);
    for (int i = 0; i < E.rows(); i++) {
        E(i, 0) = edgeSet[i].first;
        E(i, 1) = edgeSet[i].second;
    }
}

double longestEdge(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E) {
    double distance = 0;
    for (int e = 0; e < E.rows(); e++) {
        int end0 = E(e, 0);
        int end1 = E(e, 1);

        Eigen::Vector3d point0 = V.row(end0).transpose();
        Eigen::Vector3d point1 = V.row(end1).transpose();

        double localDis = (point0 - point1).norm();
        if (localDis > distance) {
            distance = localDis;
        }
    }

    return distance;
}