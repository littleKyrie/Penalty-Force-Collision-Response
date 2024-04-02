//
// Created by Kyrie Zhang on 2023/12/4.
//

#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <EigenTypes.h>
#include <Eigen/Dense>

// AABB bounding box in 3D space
template <typename T>
class BoundingBox {
    // data
public:
    Eigen::Matrix<T, 3, 1> lowerCorner;

    Eigen::Matrix<T, 3, 1> upperCorner;

    // constructors
public:
    // AABB for the entire space
    BoundingBox();

    // AABB for an edge
    BoundingBox(const Eigen::Matrix<T, 3, 1>& point1, const Eigen::Matrix<T, 3, 1>& point2);

    // AABB for a triangle
    BoundingBox(const Eigen::Matrix<T, 3, 1>& point1, const Eigen::Matrix<T, 3, 1>& point2, const Eigen::Matrix<T, 3, 1>& point3);

    // AABb for a tetrahedron
    BoundingBox(
            const Eigen::Matrix<T, 3, 1>& point1,
            const Eigen::Matrix<T, 3, 1>& point2,
            const Eigen::Matrix<T, 3, 1>& point3,
            const Eigen::Matrix<T, 3, 1>& point4);

    // AABB for an entire model
    BoundingBox(const Eigen::Matrix<T, Eigen::Dynamic, 3>& modelVertex);

    // copy constructor
    BoundingBox(const BoundingBox& other);

    // functions
public:
    // x-axis
    T width() const;

    // z-axis
    T height() const;

    // y-axis
    T depth() const;

    // the length of the appointed axis
    T length(size_t axis) const;

    T diagonalLength() const;

    T diagonalLengthSquared() const;

    Eigen::Matrix<T, 3, 1> midPoint() const;

    bool contains(const Eigen::Matrix<T, 3, 1>& point) const;

    bool overlaps(const BoundingBox& other) const;

    // reset the box to initial state(min=infinite, max=infinite)
    void reset();

    void merge(const Eigen::Matrix<T, 3, 1>& point);

    void merge(const BoundingBox& other);

    void expand(T delta);

    void scaleBox(T scaler);

    T scaleLongestEdge(T targetLength);
};

typedef BoundingBox<double> BoundingBoxD;

template <typename T>
BoundingBox<T>::BoundingBox() {
    reset();
};

template <typename T>
BoundingBox<T>::BoundingBox(const Eigen::Matrix<T, 3, 1> &point1, const Eigen::Matrix<T, 3, 1> &point2) {
    for (size_t i = 0; i < 3; i++) {
        if (point1[i] < point2[i]) {
            lowerCorner[i] = point1[i];
            upperCorner[i] = point2[i];
        } else {
            lowerCorner[i] = point2[i];
            upperCorner[i] = point1[i];
        }
    }
}

template <typename T>
BoundingBox<T>::BoundingBox(
        const Eigen::Matrix<T, 3, 1> &point1,
        const Eigen::Matrix<T, 3, 1> &point2,
        const Eigen::Matrix<T, 3, 1> &point3) {
    for (size_t i = 0; i < 3; i++) {
        T min = std::min(point1[i], point2[i]);
        lowerCorner[i] = std::min(point3[i], min);
        T max = std::max(point1[i], point2[i]);
        upperCorner[i] = std::max(point3[i], max);
    }
}

template <typename T>
BoundingBox<T>::BoundingBox(
        const Eigen::Matrix<T, 3, 1> &point1,
        const Eigen::Matrix<T, 3, 1> &point2,
        const Eigen::Matrix<T, 3, 1> &point3,
        const Eigen::Matrix<T, 3, 1> &point4) {
    for (size_t i = 0; i < 3; i++) {
        lowerCorner[i] = std::min(std::min(point1[i], point2[i]), std::min(point3[i], point4[i]));
        upperCorner[i] = std::max(std::max(point1[i], point2[i]), std::max(point3[i], point4[i]));
    }
}

template <typename T>
BoundingBox<T>::BoundingBox(const Eigen::Matrix<T, Eigen::Dynamic, 3> &modelVertex) {
    // for x
    T minX = std::numeric_limits<T>::max();
    T maxX = std::numeric_limits<T>::min();
    // for y
    T minY = std::numeric_limits<T>::max();
    T maxY = std::numeric_limits<T>::min();
    // for z
    T minZ = std::numeric_limits<T>::max();
    T maxZ = std::numeric_limits<T>::min();

    for (int i = 0; i < modelVertex.rows(); i++) {
        minX = std::min(modelVertex(i, 0), minX);
        maxX = std::max(modelVertex(i, 0), maxX);
        minY = std::min(modelVertex(i, 1), minY);
        maxY = std::max(modelVertex(i, 1), maxY);
        minZ = std::min(modelVertex(i, 2), minZ);
        maxZ = std::max(modelVertex(i, 2), maxZ);
    }

    lowerCorner = Eigen::Matrix<T, 3, 1>(minX, minY, minZ);
    upperCorner = Eigen::Matrix<T, 3, 1>(maxX, maxY, maxZ);
}

template <typename T>
BoundingBox<T>::BoundingBox(const BoundingBox<T> &other) {
    lowerCorner = other.lowerCorner;
    upperCorner = other.upperCorner;
}

template <typename T>
T BoundingBox<T>::width() const {
    return (upperCorner[0] - lowerCorner[0]);
}

template <typename T>
T BoundingBox<T>::depth() const {
    return (upperCorner[1] - lowerCorner[1]);
}

template <typename T>
T BoundingBox<T>::height() const {
    return (upperCorner[2] - lowerCorner[2]);
}

template <typename T>
T BoundingBox<T>::length(size_t axis) const {
    return (upperCorner[axis] - lowerCorner[axis]);
}

template <typename T>
T BoundingBox<T>::diagonalLength() const {
    Eigen::Matrix<T, 3, 1> diagonal = upperCorner - lowerCorner;
    return diagonal.norm();
}

template <typename T>
T BoundingBox<T>::diagonalLengthSquared() const {
    Eigen::Matrix<T, 3, 1> diagonal = upperCorner - lowerCorner;
    return diagonal.squaredNorm();
}

template <typename T>
Eigen::Matrix<T, 3, 1> BoundingBox<T>::midPoint() const {
    return Eigen::Matrix<T, 3, 1>((upperCorner[0] - lowerCorner[0]) / 2., (upperCorner[1] - lowerCorner[1]) / 2., (upperCorner[2] - lowerCorner[2]) / 2.);
}

template <typename T>
bool BoundingBox<T>::contains(const Eigen::Matrix<T, 3, 1> &point) const {
    if (point[0] < lowerCorner[0] || point[0] > upperCorner[0]) {
        return false;
    }
    if (point[1] < lowerCorner[1] || point[1] > upperCorner[1]) {
        return false;
    }
    if (point[2] < lowerCorner[2] || point[2] > upperCorner[2]) {
        return false;
    }

    return true;
}

template <typename T>
bool BoundingBox<T>::overlaps(const BoundingBox<T> &other) const {
    if (lowerCorner[0] > other.upperCorner[0] || upperCorner[0] < other.lowerCorner[0]) {
        return false;
    }
    if (lowerCorner[1] > other.upperCorner[1] || upperCorner[1] < other.lowerCorner[1]) {
        return false;
    }
    if (lowerCorner[2] > other.upperCorner[2] || upperCorner[2] < other.lowerCorner[2]) {
        return false;
    }

    return true;
}

template <typename T>
void BoundingBox<T>::reset() {
    lowerCorner[0] = std::numeric_limits<T>::min();
    lowerCorner[1] = std::numeric_limits<T>::min();
    lowerCorner[2] = std::numeric_limits<T>::min();
    upperCorner[0] = std::numeric_limits<T>::max();
    upperCorner[1] = std::numeric_limits<T>::max();
    upperCorner[2] = std::numeric_limits<T>::max();
}

template <typename T>
void BoundingBox<T>::merge(const Eigen::Matrix<T, 3, 1> &point) {
    lowerCorner[0] = std::min(lowerCorner[0], point[0]);
    lowerCorner[1] = std::min(lowerCorner[1], point[1]);
    lowerCorner[2] = std::min(lowerCorner[2], point[2]);
    upperCorner[0] = std::max(upperCorner[0], point[0]);
    upperCorner[1] = std::max(upperCorner[1], point[1]);
    upperCorner[2] = std::max(upperCorner[2], point[2]);
}

template <typename T>
void BoundingBox<T>::merge(const BoundingBox<T> &other) {
    lowerCorner[0] = std::min(lowerCorner[0], other.lowerCorner[0]);
    lowerCorner[1] = std::min(lowerCorner[1], other.lowerCorner[1]);
    lowerCorner[2] = std::min(lowerCorner[2], other.lowerCorner[2]);
    upperCorner[0] = std::max(upperCorner[0], other.upperCorner[0]);
    upperCorner[1] = std::max(upperCorner[1], other.upperCorner[1]);
    upperCorner[2] = std::max(upperCorner[2], other.upperCorner[2]);
}

template <typename T>
void BoundingBox<T>::expand(T delta) {
    for (size_t i = 0; i < 3; i++) {
        lowerCorner[i] -= delta;
        upperCorner[i] += delta;
    }
}

template <typename T>
void BoundingBox<T>::scaleBox(T scaler) {

}

template <typename T>
T BoundingBox<T>::scaleLongestEdge(T targetLength) {
    T maxLength = 0;
    maxLength = std::max(maxLength, height());
    maxLength = std::max(maxLength, width());
    maxLength = std::max(maxLength, depth());

    return static_cast<T>(targetLength / maxLength);
}

#endif //BOUNDING_BOX_H
