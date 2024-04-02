//
// Created by Kyrie Zhang on 2023/12/7.
//

#include <soft_body.h>
#include <igl/volume.h>

SoftBody::SoftBody() = default;

SoftBody::SoftBody(
        const Eigen::MatrixXd& Vertex,
        const Eigen::MatrixXi& Face,
        const Eigen::MatrixXi& Volume,
        const Eigen::MatrixXi& Edge){
    V.resize(Vertex.rows(), Vertex.cols());
    V = Vertex;
    F.resize(Face.rows(), Face.cols());
    F = Face;
    T.resize(Volume.rows(), Volume.cols());
    T = Volume;
    E.resize(Edge.rows(), Edge.cols());
    E = Edge;
    igl::volume(V, T, volumes);
}

SoftBody::SoftBody(const Eigen::MatrixXd &Vertex, const Eigen::MatrixXi &Face, const Eigen::MatrixXi &Volume) {
    V.resize(Vertex.rows(), Vertex.cols());
    V = Vertex;
    F.resize(Face.rows(), Face.cols());
    F = Face;
    T.resize(Volume.rows(), Volume.cols());
    T = Volume;
    igl::volume(V, T, volumes);
}

SoftBody::SoftBody(const SoftBody &other) {
    V.resize(other.V.rows(), other.V.cols());
    V = other.V;
    F.resize(other.F.rows(), other.F.cols());
    F = other.F;
    T.resize(other.T.rows(), other.T.cols());
    T = other.T;
    E.resize(other.E.rows(), other.E.cols());
    E = other.E;
    volumes.resize(other.volumes.size());
    volumes = other.volumes;
}

SoftBody& SoftBody::operator=(const SoftBody &other) {
    V.resize(other.V.rows(), other.V.cols());
    V = other.V;
    F.resize(other.F.rows(), other.F.cols());
    F = other.F;
    T.resize(other.T.rows(), other.T.cols());
    T = other.T;
    E.resize(other.E.rows(), other.E.cols());
    E = other.E;

    return *this;
}

void SoftBody::translate(Eigen::Ref<const Eigen::Vector3d> translation) {
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            V(i, j) += translation[j];
        }
    }
}

void SoftBody::initSelfState(Eigen::VectorXd &q, Eigen::VectorXd& qdot) const {
    q.resize(V.rows()*V.cols());
    qdot.resize(V.rows()*V.cols());

    Eigen::MatrixXd Vt = V.transpose();
    q = Eigen::Map<Eigen::VectorXd>(Vt.data(), Vt.rows()*Vt.cols());
    qdot.setZero();
}

void SoftBody::getEdgeSetForVolume() {
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

void SoftBody::getEdgeSetForSurface() {
    std::vector<std::pair<int, int>> tmp_edge;
    for (int t = 0; t < F.rows(); t++) {
        Eigen::RowVector3i tet = F.row(t);
        int p0 = tet[0];
        int p1 = tet[1];
        int p2 = tet[2];

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
        if (p1 < p2) {
            tmp_edge.emplace_back(p1, p2);
        } else {
            tmp_edge.emplace_back(p2, p1);
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

bool SoftBody::isTriangleEdge(
        Eigen::Ref<const Eigen::RowVector2i> edge,
        Eigen::Ref<const Eigen::RowVector3i> triangle) {
    int e0 = edge[0];
    int e1 = edge[1];
    int t0 = triangle[0];
    int t1 = triangle[1];
    int t2 = triangle[2];

    // (t0, t1) & (e0, e1)
    if ((e0 == t0 && e1 == t1) || (e0 == t1 && e1 == t0)) {
        return true;
    }
    // (t0, t2) & (e0, e1)
    if ((e0 == t0 && e1 == t2) || (e0 == t2 && e1 == t0)) {
        return true;
    }
    // (t1, t2) & (e0, e1)
    if ((e0 == t1 && e1 == t2) || (e0 == t2 && e1 == t1)) {
        return true;
    }

    return false;
}

bool SoftBody::isTetrahedralEdge(
        Eigen::Ref<const Eigen::RowVector2i> edge,
        Eigen::Ref<const Eigen::RowVector4i> tet) {
    Eigen::RowVector3i tri0(tet[0], tet[1], tet[2]);
    Eigen::RowVector3i tri1(tet[0], tet[1], tet[3]);
    Eigen::RowVector3i tri2(tet[0], tet[2], tet[3]);
    Eigen::RowVector3i tri3(tet[1], tet[2], tet[3]);

    if (isTriangleEdge(edge, tri0)) {
        return true;
    }
    if (isTriangleEdge(edge, tri1)) {
        return true;
    }
    if (isTriangleEdge(edge, tri2)) {
        return true;
    }
    if (isTriangleEdge(edge, tri3)) {
        return true;
    }

    return false;
}

bool SoftBody::isEdgeConnectTriangle(
        Eigen::Ref<const Eigen::RowVector2i> edge,
        Eigen::Ref<const Eigen::RowVector3i> tri) {
    int e0 = edge[0];
    int e1 = edge[1];

    int t0 = tri[0];
    int t1 = tri[1];
    int t2 = tri[2];

    if (e0 == t0 || e0 == t1 || e0 == t2 || e1 == t0 || e1 == t1 || e1 == t2) {
        return true;
    }

    return false;
}