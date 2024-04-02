//
// Created by Kyrie Zhang on 2023/12/4.
//

#include <spatial_hash.h>

SpatialHash::SpatialHash() {
    initialize();
}

SpatialHash::SpatialHash(int numX, int numY, int numZ, double gridSpacing) {
    resize(numX, numY, numZ);
    reset(gridSpacing);
}

SpatialHash::SpatialHash(int numX, int numY, int numZ) {
    resize(numX, numY, numZ);
}

SpatialHash::SpatialHash(double gridSpacing) {
    initialize();
    reset(gridSpacing);
}

SpatialHash::SpatialHash(const SpatialHash &other) {
    _x = other._x;
    _y = other._y;
    _z = other._z;
    _gridSpacing = other._gridSpacing;
    _cells = other._cells;
}

void SpatialHash::clear() {
    for (auto &_cell : _cells) {
        _cell.clear();
    }
}

void SpatialHash::initialize() {
    _cells.resize(_x*_y*_z);
    clear();
}

void SpatialHash::resize(int numX, int numY, int numZ) {
    resetX(numX);
    resetY(numY);
    resetZ(numZ);
    initialize();
}

void SpatialHash::reset(double gridSpacing) {
    _gridSpacing = gridSpacing;
}

void SpatialHash::resetX(int numX) {
    assert(numX > 0);
    _x = numX;
}

void SpatialHash::resetY(int numY) {
    assert(numY > 0);
    _y = numY;
}

void SpatialHash::resetZ(int numZ) {
    assert(numZ > 0);
    _z = numZ;
}

int SpatialHash::getX() const {
    return _x;
}

int SpatialHash::getY() const {
    return _y;
}

int SpatialHash::getZ() const {
    return _z;
}

double SpatialHash::getGridSize() const {
    return _gridSpacing;
}

std::vector<std::pair<int, int>>& SpatialHash::operator()(int i, int j, int k) {
    // (i, j, k) represents the index of global spatial grid
    Eigen::Vector3i cellIndex(i, j, k);

    // use wrapped index to transfer the index to the hash space
    Eigen::Vector3i wrappedIndex = getHashIndexFromCellIndex(cellIndex);

    size_t key = wrappedIndex[1] * (_x * _z) + wrappedIndex[2] * _x + wrappedIndex[0];
    return _cells[key];
}

const std::vector<std::pair<int, int>>& SpatialHash::operator()(int i, int j, int k) const {
    // (i, j, k) represents the index of global spatial grid
    Eigen::Vector3i cellIndex(i, j, k);

    // use wrapped index to transfer the index to the hash space
    Eigen::Vector3i wrappedIndex = getHashIndexFromCellIndex(cellIndex);

    size_t key = wrappedIndex[1] * (_x * _z) + wrappedIndex[2] * _x + wrappedIndex[0];
    return _cells[key];
}

std::vector<std::pair<int, int>>& SpatialHash::operator[](size_t key) {
    assert(key < _cells.size());

    return _cells[key];
}

const std::vector<std::pair<int, int>>& SpatialHash::operator[](size_t key) const {
    assert(key < _cells.size());

    return _cells[key];
}

SpatialHash& SpatialHash::operator=(const SpatialHash &table) {
    _cells = table._cells;
    _x = table._x;
    _y = table._y;
    _z = table._z;
    _gridSpacing = table._gridSpacing;

    return *this;
}

Eigen::Vector3i SpatialHash::keyToIndex(size_t key) const {
    assert(key < _cells.size());

    // wrapped index
    Eigen::Vector3i wrappedIndex;
    wrappedIndex[0] = key % _x;
    wrappedIndex[1] = (key / _x) / _z;
    wrappedIndex[2] = (key / _x) % _z;

    return wrappedIndex;
}

Eigen::Vector3i SpatialHash::getCellIndexFromNodePosition(const Eigen::Vector3d& node) const {
    Eigen::Vector3i cellIndex;
    cellIndex[0] = int(std::floor(node.x() / _gridSpacing));
    cellIndex[1] = int(std::floor(node.y() / _gridSpacing));
    cellIndex[2] = int(std::floor(node.z() / _gridSpacing));

    return cellIndex;
}

Eigen::Vector3i SpatialHash::getHashIndexFromCellIndex(const Eigen::Vector3i &cellIndex) const {
    Eigen::Vector3i wrappedIndex = cellIndex;
    wrappedIndex[0] = cellIndex[0] % _x;
    wrappedIndex[1] = cellIndex[1] % _y;
    wrappedIndex[2] = cellIndex[2] % _z;

    if (wrappedIndex[0] < 0) {
        wrappedIndex[0] += _x;
    }
    if (wrappedIndex[1] < 0) {
        wrappedIndex[1] += _y;
    }
    if (wrappedIndex[2] < 0) {
        wrappedIndex[2] += _z;
    }

    return wrappedIndex;
}

size_t SpatialHash::getHashKeyFromNodePosition(const Eigen::Vector3d& node) const {
    Eigen::Vector3i cellIndex = getCellIndexFromNodePosition(node);

    // map key from global spatial grid to hash space
    Eigen::Vector3i wrappedIndex = getHashIndexFromCellIndex(cellIndex);

    size_t key = wrappedIndex[1] * _x * _z + wrappedIndex[2] * _x + wrappedIndex[0];
    return key;
}

void SpatialHash::build(const Eigen::VectorXd& positions, int id) {
    size_t numVertex = positions.size() / 3;
    for (int i = 0; i < numVertex; i++) {
        Eigen::Vector3d point = positions.segment(3*i, 3);
        size_t key = getHashKeyFromNodePosition(point);
        _cells[key].emplace_back(i, id);
    }
}

void SpatialHash::build(
        const Eigen::MatrixXi& edgeSet,
        const Eigen::VectorXd& positions,
        int id) {
    for (int i = 0; i < edgeSet.rows(); i++) {
        Eigen::RowVector2i edge = edgeSet.row(i);
        int p0 = edge[0];
        int p1 = edge[1];

        Eigen::Vector3d node0 = positions.segment(3*p0, 3);
        Eigen::Vector3d node1 = positions.segment(3*p1, 3);
        BoundingBoxD AABB(node0, node1);
        std::vector<size_t> keys;
        coveredCells(keys, AABB);
        for (const auto &key : keys) {
            _cells[key].emplace_back(i, id);
        }
    }
}

void SpatialHash::coveredCells(std::vector<size_t>& keyList, const BoundingBoxD &AABB) const {
    // get the range of cell indices
    Eigen::Vector3i lowerCornerIndex = getCellIndexFromNodePosition(AABB.lowerCorner);
    Eigen::Vector3i upperCornerIndex = getCellIndexFromNodePosition(AABB.upperCorner);

    // get all covered cell indices
    for (int i = lowerCornerIndex[0]; i <= upperCornerIndex[0]; i++)
        for (int k = lowerCornerIndex[2]; k <= upperCornerIndex[2]; k++)
            for (int j = lowerCornerIndex[1]; j <= upperCornerIndex[1]; j++) {
                // wrap the cell index
                Eigen::Vector3i wrappedIndex = getHashIndexFromCellIndex(Eigen::Vector3i(i, j, k));
                size_t key = wrappedIndex[1] * _x * _z + wrappedIndex[2] * _x + wrappedIndex[0];
                keyList.push_back(key);
            }
}

void SpatialHash::boundaryKeys(Boundary &boundary) const {
    // bottom & top
    if (boundary.normal[1] == 1 || boundary.normal[1] == -1) {
        double y = boundary.point[1];
        int j = int(std::floor(y / _gridSpacing));
        for (int i = 0; i < _x; i++)
            for (int k = 0; k < _z; k++) {
                Eigen::Vector3i wrappedIndex = getHashIndexFromCellIndex(Eigen::Vector3i(i, j, k));
                size_t key = wrappedIndex[1] * (_x * _z) + wrappedIndex[2] * _x + wrappedIndex[0];
                boundary.hashKeys.push_back(key);
            }

        return;
    }

    // left & right
    if (boundary.normal[0] == 1 || boundary.normal[0] == -1) {
        double x = boundary.point[0];
        int i = int(std::floor(x / _gridSpacing));
        for (int k = 0; k < _z; k++)
            for (int j = 0; j < _y; j++) {
                Eigen::Vector3i wrappedIndex = getHashIndexFromCellIndex(Eigen::Vector3i(i, j, k));
                size_t key = wrappedIndex[1] * (_x * _z) + wrappedIndex[2] * _x + wrappedIndex[0];
                boundary.hashKeys.push_back(key);
            }

        return;
    }

    // forward & backward
    if (boundary.normal[2] == -1 || boundary.normal[2] == 1) {
        double z = boundary.point[2];
        int k = int(std::floor(z / _gridSpacing));
        for (int i = 0; i < _x; i++)
            for (int j = 0; j < _y; j++) {
                Eigen::Vector3i wrappedIndex = getHashIndexFromCellIndex(Eigen::Vector3i(i, j, k));
                size_t key = wrappedIndex[1] * (_x * _z) + wrappedIndex[2] * _x + wrappedIndex[0];
                boundary.hashKeys.push_back(key);
            }

        return;
    }

}

// for test
void SpatialHash::visualHashTable() const {
    for (int key = 0; key < _cells.size(); key++) {
        if (!_cells[key].empty()) {
            Eigen::Vector3i index = keyToIndex(key);
            std::cout << "( " << index.x() << ", " << index.y() << ", " << index.z() << " ): "<< std::endl;
            for (const auto &pair : _cells[key]) {
                std::cout << "at obj " << pair.second << " has vertex " << pair.first << ", ";
            }
            std::cout << std::endl;
        }
    }
}