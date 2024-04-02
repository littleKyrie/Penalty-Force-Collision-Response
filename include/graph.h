//
// Created by Kyrie Zhang on 2023/11/16.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <thread>
#include <vector>
#include <memory>

#include <oneapi/tbb.h>
#include <SparseStorage.h>

#include <EigenTypes.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

/*
 * create a simple undirected graph for parallel Gauss Seidel
 */
class Vertex {
private:
    int number;
    int color;
    bool isVisited;
    bool isInQueue;

    // constructors
public:
    Vertex();

    explicit Vertex(int index);

    Vertex(int index, size_t colour);

    Vertex(const Vertex& other);

    // setter and getters
public:
    void setColor(size_t colour);

    int getNumber() const;

    int getColor() const;

    bool visited() const;

    void visitVertex();

    bool inQueue() const;

    void queueVertex();

    // functions
public:
    bool operator==(const Vertex& other) const;
};

typedef std::shared_ptr<Vertex> VertexPtr;

class Edge {
    // data
private:
    Vertex vertex1;
    Vertex vertex2;

    // constructors
public:
    Edge();

    Edge(const Vertex& v1, const Vertex& v2);

    // getters
public:
    Vertex x() const;

    Vertex y() const;

    // functions
public:
    bool operator==(const Edge& other) const;
};

class Palette {
    //data
private:
    std::vector<size_t> pigment;
    std::vector<bool> isUsed;

    // constructor
public:
    Palette();

    Palette(const std::vector<size_t>& colors);

    // getter & setter
public:
    size_t getNumberOfPigment() const;

    void addColor();

    // functions
public:
    bool availablePigment() const;

    void fillPigments();

    void coloration(Vertex& v);

    void colorByNeighbors(VertexPtr& v, const std::vector<VertexPtr>& neighbors);
};

class Graph {
    // data
private:
    std::vector<Vertex> V;
    std::vector<Edge> E;
    CCSi adjacencyMatrix;

    // constructor
public:
    Graph();

    // resize
public:
    void resizeV(size_t numV);

    void resizeE(size_t numE);

    void clear();

    size_t sizeV() const;

    size_t sizeE() const;

    // storage
public:
    // concurrent implementation
    void buildGraph(const std::vector<Vertex>& vertexSet,
                    const std::vector<Edge>& edgeSet);

    // functions
public:
    bool adjacent(const Vertex& v1,
                  const Vertex& v2) const;

    void neighbors(
            std::vector<VertexPtr>& neighbor,
            const VertexPtr& v) const;

    size_t degreeOfGraph() const;

    // color the vertex
    size_t broadFirstSearch();

    int visitVertexColor(size_t index) const;

    int visitVertexNumber(size_t index) const;

    // concurrent
    void parallelColorGraph();

};

#endif //GRAPH_H
