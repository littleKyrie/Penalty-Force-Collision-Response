//
// Created by Kyrie Zhang on 2023/12/4.
//

#ifndef SPATIAL_HASH_H
#define SPATIAL_HASH_H

#include <vector>
#include <memory>
#include <bounding_box.h>
#include <boundary.h>
#include <iostream>

// from the paper <Optimized Spatial Hashing for Collision Detection of Deformable>
class SpatialHash {
    //data
private:
    // hash table
    std::vector<std::vector<std::pair<int, int>>> _cells; // <number of primitive, id of object>
    // the coordinate system is the same with the OpenGL
    int _x = 1; // the number of cells in x-axis
    int _y = 1; // the number of cells in y-axis
    int _z = 1; // the number of cells in z-axis
    double _gridSpacing = 1.0; // the size of a single cell

    // constructors
public:
    // use default resolution, grid size and origin to initialize the spatial grid
    SpatialHash();

    // reset the resolution and grid size to initialize the grid
    SpatialHash(int numX, int numY, int numZ, double gridSpacing);

    // reset the resolution
    SpatialHash(int numX, int numY, int numZ);

    // reset the grid size
    explicit SpatialHash(double gridSpacing);

    SpatialHash(const SpatialHash& other);

    // functions
public:
    // make the hash table no elements in each key
    void clear();

    // initialize the size of hash table according to current space resolution and clear elements in each cell
    void initialize();

    // regulate the number of cells in each axis (regulate the size of hash table) and clear elements
    void resize(int numX, int numY, int numZ);

    // reset the size of a single cell
    void reset(double gridSpacing);

    std::vector<std::pair<int, int>>& operator()(int i, int j, int k);

    const std::vector<std::pair<int, int>>& operator()(int i, int j, int k) const;

    std::vector<std::pair<int, int>>& operator[](size_t key);

    const std::vector<std::pair<int, int>>& operator[](size_t key) const;

    SpatialHash& operator=(const SpatialHash& table);

    Eigen::Vector3i keyToIndex(size_t key) const;

    // only for vertices of a single object
    void build(const Eigen::VectorXd& positions,int id);

    // only for edge of a single object
    void build(const Eigen::MatrixXi& edgeSet, const Eigen::VectorXd& positions, int id);

    Eigen::Vector3i getCellIndexFromNodePosition(const Eigen::Vector3d& node) const;

    Eigen::Vector3i getHashIndexFromCellIndex(const Eigen::Vector3i& cellIndex) const;

    size_t getHashKeyFromNodePosition(const Eigen::Vector3d& node) const;

    void coveredCells(
            std::vector<size_t>& keyList,
            const BoundingBoxD& AABB) const;

    // I use planes which are parallel with xoy, xoz or yoz as the boundary and use boundary to create the grid,
    // so the boundary covered cells are easy to define
    void boundaryKeys(Boundary& boundary) const;

    // getters
public:
    int getX() const;

    int getY() const;

    int getZ() const;

    double getGridSize() const;

    // functions
protected:
    void resetX(int numX);

    void resetY(int numY);

    void resetZ(int numZ);

    // for test
public:
    void visualHashTable() const;
};

#endif //SPATIAL_HASH_H