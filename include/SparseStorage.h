//
// Created by Kyrie Zhang on 2023/11/14.
//

#ifndef SPARSESTORAGE_H
#define SPARSESTORAGE_H

#include <iostream>
#include <memory>
#include <thread>
#include <vector>

#include <oneapi/tbb.h>

template <typename dataType>
class CompressedColumnStorage {
    // typedef
public:
    typedef std::tuple<size_t, size_t, dataType> triplet;
    typedef tbb::concurrent_vector<triplet> tripletList;

    // data
private:
    std::vector<dataType> values;
    std::vector<size_t> innerIndices;// row index of any none zero element
    std::vector<size_t> outerStarts; // the index of the first none zero element of each column in previous two arrays
    size_t rowSize;
    size_t colSize;

    // constructors
public:
    CompressedColumnStorage();

    CompressedColumnStorage(size_t r, size_t c);

public:
    void resize(size_t newRow, size_t newCol);

    size_t rows() const;

    size_t cols() const;

    void clear();

public:
    bool insertValue(size_t i, size_t j, dataType value);

    dataType coeff(size_t i, size_t j) const;

    double* coeffRef(size_t i, size_t j);

    // concurrent constructing a complete sparse structure by a column list
    // which is a vector with size of the number of columns
    // and has elements consisted of a triplet list
    // which include some triplet with the index of matrix element in current column
    // and relevant value
    bool setFromColumnList(const std::vector<tripletList>& columnList);

    bool setFromColumnClusters(const std::vector<tbb::concurrent_vector<dataType>>& columnClusters);
};

typedef CompressedColumnStorage<double> CCSd;
typedef CompressedColumnStorage<float> CCSf;
typedef CompressedColumnStorage<int> CCSi;

template <typename dataType>
CompressedColumnStorage<dataType>::CompressedColumnStorage() {
    rowSize = 0;
    colSize = 0;
};

template <typename dataType>
CompressedColumnStorage<dataType>::CompressedColumnStorage(size_t newRow, size_t newCol) {
    resize(newRow, newCol);
}

template <typename dataType>
void CompressedColumnStorage<dataType>::resize(size_t newRow, size_t newCol) {
    if (newRow == 0 || newCol == 0) {
        std::cout << "new size of row and col should > 0" << std::endl;
        exit(0);
    }

    rowSize = newRow;
    colSize = newCol;
}

template <typename dataType>
size_t CompressedColumnStorage<dataType>::rows() const {
    return rowSize;
}

template <typename dataType>
size_t CompressedColumnStorage<dataType>::cols() const {
    return colSize;
}

template <typename dataType>
void CompressedColumnStorage<dataType>::clear() {
    rowSize = 0;
    colSize = 0;

    values.clear();
    innerIndices.clear();
    outerStarts.clear();
}

template <typename dataType>
bool CompressedColumnStorage<dataType>::insertValue(size_t i, size_t j, dataType value) {
    if (i >= rowSize || j >= colSize) {
        std::cout << "the access is out of bounds" << std::endl;
        return false;
    }

    size_t start = outerStarts[j];
    size_t end = outerStarts[j+1];
    // binary search
    int low = start;
    int high = end - 1;
    int mid = 0;
    while (low <= high) {
        mid = (low + high) / 2;
        if (innerIndices[mid] == i) {
            std::cout << "the value has existed, you should use coeffRef to change the value" << std::endl;
            return false;
        } else if (innerIndices[mid] > i) {
            high = mid - 1;
        } else if (innerIndices[mid] < i) {
            low = mid + 1;
        }
    }

    if (value == 0) {
        return true;
    }

    values.insert(values.begin()+mid, value);
    innerIndices.insert(innerIndices.begin()+mid, i);
    for (size_t k = j + 1; k < outerStarts.size(); k++) {
        ++outerStarts[k];
    }

    // after insert for test
    std::cout << "after insert" << std::endl;
    std::cout << "values" << std::endl;
    for (size_t i = 0; i < values.size(); i++) {
        std::cout << values[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "inner indices" << std::endl;
    for (size_t i = 0; i < innerIndices.size(); i++) {
        std::cout << innerIndices[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "outerStarts" << std::endl;
    for (size_t i = 0; i < outerStarts.size(); i++) {
        std::cout << outerStarts[i] << ", ";
    }
    std::cout << std::endl;

    return true;
}

template <typename dataType>
dataType CompressedColumnStorage<dataType>::coeff(size_t i, size_t j) const {
    if (i >= rowSize || j >= colSize) {
        std::cout << "the access is out of bounds" << std::endl;
        exit(0);
    }

    size_t start = outerStarts[j];
    size_t end = outerStarts[j+1];
    // binary search
    int low = start;
    int high = end - 1;
    int mid = 0;
    while (low <= high) {
        mid = (low + high) / 2;
        if (innerIndices[mid] == i) {
            dataType result = values[mid];
            return result;
        } else if (innerIndices[mid] > i) {
            high = mid - 1;
        } else if (innerIndices[mid] < i) {
            low = mid + 1;
        }
    }
    return 0;
}

template <typename dataType>
double* CompressedColumnStorage<dataType>::coeffRef(size_t i, size_t j) {
    if (i >= rowSize || j >= colSize) {
        std::cout << "the access is out of bounds" << std::endl;
        exit(0);
    }

    size_t start = outerStarts[j];
    size_t end = outerStarts[j+1];
    // binary search
    int low = start;
    int high = end - 1;
    int mid = 0;
    while (low <= high) {
        mid = (low + high) / 2;
        if (innerIndices[mid] == i) {
            return values[mid];
        } else if (innerIndices[mid] > i) {
            high = mid - 1;
        } else if (innerIndices[mid] < i) {
            low = mid + 1;
        }
    }

    std::cout << "the value doesn't exist, you should insert a none zero value first" << std::endl;
    return nullptr;
}

template <typename dataType>
bool CompressedColumnStorage<dataType>::setFromColumnList(const std::vector<tripletList> &columnList) {
    // in this function, I just consider the special matrix
    // which means no column only have 0 as its element
    size_t n = columnList.size();
    if (n < rowSize) {
        std::cout << "the matrix has 0 column" << std::endl;
        return false;
    }
    outerStarts.resize(n+1);
    outerStarts[0] = 0;
    tbb::parallel_for(tbb::blocked_range<size_t>(1, outerStarts.size()),
                      [&](const tbb::blocked_range<size_t>& r){
                          for (size_t i = r.begin(); i < r.end(); i++) {
                              size_t sum = 0;
                              for (size_t j = 0; j < i; j++) {
                                  sum += columnList[j].size();
                              }
                              outerStarts[i] = sum;
                          }
                      });

    size_t nnz = outerStarts[n];
    values.resize(nnz);
    innerIndices.resize(nnz);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
                      [&](const tbb::blocked_range<size_t>& r) {
                          for (size_t i = r.begin(); i < r.end(); i++) {
                              tripletList list = columnList[i];
                              // sort
                              std::sort(list.begin(), list.end(),
                                        [](const triplet& tuple1, const triplet & tuple2) {
                                            return std::get<0>(tuple1) < std::get<0>(tuple2);
                                        });

                              size_t start = outerStarts[i];
                              size_t end = outerStarts[i+1];
                              for (size_t j = 0; j < end - start; j++) {
                                  values[start+j] = std::get<2>(list[j]);
                                  innerIndices[start+j] = std::get<0>(list[j]);
                              }
                          }
                      });

    return true;
}

template <typename dataType>
bool CompressedColumnStorage<dataType>::setFromColumnClusters(
        const std::vector<tbb::concurrent_vector<dataType>> &columnClusters) {
    std::vector<size_t> nnzOfEachCol(colSize);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, colSize),
            [&](const tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i < r.end(); i++) {
            size_t start = i * rowSize;
            size_t end = start + rowSize;
            size_t sum = 0;
            for (size_t j = start; j < end; j++) {
                if (columnClusters[j].size() != 0) {
                    ++sum;
                }
            }
            nnzOfEachCol[i] = sum;
        }
    });

    outerStarts.resize(colSize+1);
    outerStarts[0] = 0;
    size_t nnz = tbb::parallel_scan(tbb::blocked_range<size_t>(1, outerStarts.size()), 0,
            [&](const tbb::blocked_range<size_t>& r, size_t sum, bool is_final_scan) -> size_t {
            size_t tmp = sum;
            for (size_t i = r.begin(); i < r.end(); i++) {
                tmp += nnzOfEachCol[i-1];
                if (is_final_scan) {
                    outerStarts[i] = tmp;
                }
            }
            return tmp;
        },
            [](size_t left, size_t right) -> size_t {
                return left + right;
    });
    values.resize(nnz);
    innerIndices.resize(nnz);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, colSize),
            [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
            size_t start = outerStarts[i];
            size_t k = 0;
            size_t end = i * rowSize + rowSize;

            for (size_t j = i * rowSize; j < end; j++) {
                if (columnClusters[j].size() != 0) {
                    innerIndices[start+k] = j % colSize;
                    for (const auto &ele : columnClusters[j]) {
                        values[k+start] += ele;
                    }
                    k++;
                }
            }
        }
    });
}

#endif //SPARSESTORAGE_H
