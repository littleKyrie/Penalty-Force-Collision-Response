//
// Created by Kyrie Zhang on 2023/11/16.
//

#include <graph.h>
#include <queue>
#include <src/SparseStorage.cpp>

Vertex::Vertex() {
    number = -1;
    color = -1;
    isVisited = false;
    isInQueue = false;
}

Vertex::Vertex(int index) {
    number = index;
    color = -1;
    isVisited = false;
    isInQueue = false;
}

Vertex::Vertex(int index, size_t colour) {
    number = index;
    color = colour;
    isVisited = false;
    isInQueue = false;
}

Vertex::Vertex(const Vertex &other) {
    number = other.number;
    color = other.color;
    isVisited = other.isVisited;
    isInQueue = other.isInQueue;
}

int Vertex::getNumber() const {
    return number;
}

int Vertex::getColor() const {
    return color;
}

void Vertex::setColor(size_t colour) {
    color = colour;
}

bool Vertex::visited() const {
    return isVisited;
}

void Vertex::visitVertex() {
    isVisited = true;
}

bool Vertex::inQueue() const {
    return isInQueue;
}

void Vertex::queueVertex() {
    isInQueue = true;
}

bool Vertex::operator==(const Vertex &other) const {
    return number == other.number;
}

//--------------------------------
Edge::Edge() = default;

Edge::Edge(const Vertex &v1, const Vertex &v2) {
    vertex1 = v1;
    vertex2 = v2;
}

Vertex Edge::x() const {
    return vertex1;
}

Vertex Edge::y() const {
    return vertex2;
}

bool Edge::operator==(const Edge &other) const {
    if (vertex1 == other.x() && vertex2 == other.y()) {
        return true;
    } else if (vertex1 == other.y() && vertex2 == other.x()) {
        return true;
    }

    return false;
}

//--------------------------------
Palette::Palette() {
    pigment.push_back(0);
    isUsed.push_back(false);
}

Palette::Palette(const std::vector<size_t> &colors) {
    pigment.resize(colors.size());
    isUsed.resize(colors.size());
    for (size_t i = 0; i < colors.size(); i++) {
        pigment[i] = colors[i];
        isUsed[i] = false;
    }
}

size_t Palette::getNumberOfPigment() const {
    return pigment.size();
}

void Palette::addColor() {
    size_t oldColor = pigment.back();
    pigment.push_back(++oldColor);
    isUsed.push_back(false);
}

bool Palette::availablePigment() const {
    for (size_t i = 0; i < getNumberOfPigment(); i++) {
        if (!isUsed[i]) {
            return true;
        }
    }

    return false;
}

void Palette::fillPigments() {
    if (!availablePigment()) {
        addColor();
    }
}

void Palette::coloration(Vertex &v) {
    for (size_t i = 0; i < getNumberOfPigment(); i++) {
        if (!isUsed[i]) {
            v.setColor(pigment[i]);
            isUsed[i] = true;
        }
    }

    addColor();
    v.setColor(pigment.back());
    isUsed[getNumberOfPigment()-1] = true;
}

void Palette::colorByNeighbors(
        VertexPtr &v,
        const std::vector<VertexPtr> &neighbors) {
    if (neighbors.empty()) {
        v->setColor(0);
        return;
    }
    // size_t num = neighbors.size();
    for (const auto& neighbor : neighbors) {
        size_t currentColor = neighbor->getColor();
        int start = 0;
        int end = getNumberOfPigment();

        int low = start;
        int high = end - 1;
        int mid = 0;

        while (low <= high) {
            mid = (low + high) / 2;
            if (pigment[mid] == currentColor) {
                isUsed[mid] = true;
                break;
            } else if (pigment[mid] > currentColor) {
                high = mid - 1;
            } else if (pigment[mid] < currentColor) {
                low = mid + 1;
            }
        }
    }

    size_t i;
    for (i = 0; i < getNumberOfPigment(); i++) {
        if (!isUsed[i]) {
            v->setColor(pigment[i]);
            break;
        }
    }

    // for test
    /*
    std::cout << "in old palette: " << std::endl;
    for (size_t i = 0; i < getNumberOfPigment(); i++) {
        std::cout << "color " << pigment[i] << " and isUsed as " << isUsed[i] << std::endl;
    }
     */

    if (i == getNumberOfPigment()) {
        addColor();
        v->setColor(pigment.back());
    }

    // for test
    // std::cout << "in new palette: " << std::endl;
    for (size_t j = 0; j < getNumberOfPigment(); j++) {
        isUsed[j] = false;
        // std::cout << "color " << pigment[j] << " and isUsed as " << isUsed[j] << std::endl;
    }
}

//--------------------------------
Graph::Graph() {
    V.resize(0);
    E.resize(0);
    // adjacencyMatrix.resize(0, 0);
}

void Graph::resizeV(size_t numV) {
    V.resize(numV);
    adjacencyMatrix.resize(numV, numV);
}

void Graph::resizeE(size_t numE) {
    E.resize(numE);
}

void Graph::clear() {
    V.clear();
    E.clear();
    adjacencyMatrix.clear();
}

size_t Graph::sizeV() const {
    return V.size();
}

size_t Graph::sizeE() const {
    return E.size();
}

void Graph::buildGraph(
        const std::vector<Vertex> &vertexSet,
        const std::vector<Edge> &edgeSet) {
    size_t numV = vertexSet.size();
    size_t numE = edgeSet.size();

    resizeV(numV);
    resizeE(numE);
    adjacencyMatrix.resize(numV, numV);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numV),
            [&](const tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i < r.end(); i++) {
            V[i] = vertexSet[i];
        }
    });
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numE),
            [&](const tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i < r.end(); i++) {
            E[i] = edgeSet[i];
        }
    });

    std::vector<tbb::concurrent_vector<std::tuple<size_t, size_t, int>>> columnList;
    columnList.resize(adjacencyMatrix.cols());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, sizeE()),
            [&](const tbb::blocked_range<size_t>& r){
        for (size_t i = r.begin(); i < r.end(); i++) {
            size_t colIndex = E[i].y().getNumber();
            size_t rowIndex = E[i].x().getNumber();

            std::tuple<size_t, size_t, int> triplet1(rowIndex, colIndex, 1);
            columnList[colIndex].push_back(triplet1);

            // the matrix is symmetric
            std::tuple<size_t, size_t, int> triplet2(colIndex, rowIndex, 1);
            columnList[rowIndex].push_back(triplet2);
        }
    });
    adjacencyMatrix.setFromColumnList(columnList);

    // for test
    /*
    for (size_t i = 0; i < adjacencyMatrix.rows(); i++) {
        for (size_t j = 0; j < adjacencyMatrix.cols(); j++) {
            std::cout << adjacencyMatrix.coeff(i, j) << ", ";
        }
        std::cout << std::endl;
    }
     */
}

bool Graph::adjacent(
        const Vertex &v1,
        const Vertex &v2) const {
    int i = v1.getNumber();
    int j = v2.getNumber();

    if (adjacencyMatrix.coeff(i, j) == 0) {
        return false;
    } else {
        return true;
    }
}

void Graph::neighbors(
        std::vector<VertexPtr> &neighbor,
        const VertexPtr &v) const {
    int i = v->getNumber();
    for (size_t j = 0; j < adjacencyMatrix.cols(); j++) {
        if (i == j) {
            continue;
        } else if (adjacencyMatrix.coeff(i, j) == 1) {
            // if when we construct vertex set,
            // we make the vertex number = vertex index in the set
            // there can be more fast
            if (V[j].getNumber() == j) {
                VertexPtr vj = std::make_shared<Vertex>(V[j]);
                neighbor.push_back(vj);
            } else {
                std::cout << "the order in V have some problem" << std::endl;
                exit(0);
            }
        }
    }
}

size_t Graph::degreeOfGraph() const {
    /*
    tbb::parallel_for(tbb::blocked_range<size_t>(0, sizeV()), [&](const tbb::blocked_range<size_t>& r){

    });
     */
    size_t degree = 0;
    for (size_t i = 0; i < adjacencyMatrix.rows(); i++) {
        size_t sum = 0;
        for (size_t j = 0; j < adjacencyMatrix.cols(); j++) {
            if (adjacencyMatrix.coeff(i, j) == 1 && i != j) {
                sum++;
            }
        }
        if (degree < sum) {
            degree = sum;
        }
    }

    return degree;
}

size_t Graph::broadFirstSearch() {
    Palette palette;
    std::queue<VertexPtr> visitQueue;
    size_t numV = sizeV();
    if (numV == 0) {
        return 0;
    }

    // V[0].visitVertex();
    // V[0].setColor(4);
    VertexPtr v0 = std::make_shared<Vertex>(V[0]);
    visitQueue.push(v0);
    /*
    std::cout << "v0 isVisited "<< v0->visited() << std::endl;
    v0->visitVertex();
    std::cout << "change v0 isVisited " << v0->visited() << std::endl;
    std::cout << "V[0] isVisited " << V[0].visited() << std::endl;
    std::cout << "v0 color is " << v0->getColor() << std::endl;
    v0->setColor(5);
    std::cout << "v0 new color is " << v0->getColor() << std::endl;
    std::cout << "V0 color is " << V[0].getColor() << std::endl;
    V[0].setColor(v0->getColor());
    std::cout << "V0 new color is " << V[0].getColor() << std::endl;
     */

    while(!visitQueue.empty()) {
        std::vector<VertexPtr> neighbor;
        VertexPtr head = visitQueue.front();
        neighbors(neighbor, head);

        // for test
        /*
        std::cout << "the vertex " << head->getNumber() << " has neighbor:" << std::endl;
        for (const auto & vertex : neighbor) {
             std::cout << vertex->getNumber() << " ";
             std::cout << "isVisited is " << vertex->visited();
        }
        std::cout << std::endl;
        */

        palette.colorByNeighbors(head, neighbor);
        // head->visitVertex();
        // write back
        V[head->getNumber()].visitVertex();
        V[head->getNumber()].setColor(head->getColor());
        for (const auto& v : neighbor) {
            // size_t num = v->getNumber();
            if (!v->visited() && !v->inQueue()) {
                visitQueue.push(v);
                V[v->getNumber()].queueVertex();
            }
        }

        // for test
        /*
        std::cout << "after color, the vertex " << V[head->getNumber()].getNumber() << " have some new property" << std::endl;
        std::cout << "the color is " << V[head->getNumber()].getColor() << " and visited is " << V[head->getNumber()].visited() << std::endl;
        */

        visitQueue.pop();

        // for test
        /*
        std::cout << "the info of queue is: " << std::endl;
        size_t queueSize = visitQueue.size();
        std::cout << "the size of new queue is " << queueSize << std::endl;
        std::cout << "the elements in the queue are: ";
        std::vector<size_t> currentQueue;
        for (size_t a = 0; a < queueSize; a++) {
            VertexPtr ele = visitQueue.front();
            currentQueue.push_back(ele->getNumber());
            visitQueue.pop();
            visitQueue.push(ele);
        }
        for (size_t a = 0; a < currentQueue.size(); a++) {
            std::cout << currentQueue[a] << ", ";
        }
        std::cout << std::endl;
         */

        if (visitQueue.empty()) {
            for (size_t i = 0; i < numV; i++) {
                if (!V[i].visited()) {
                    VertexPtr vi = std::make_shared<Vertex>(V[i]);
                    visitQueue.push(vi);
                }
            }
        }
    }

    // test if some nodes has the same color with their neighbors
    /*
    std::vector<std::pair<int, int>> result;
    // V[1].setColor(0);
    for (size_t i = 0; i < V.size(); i++) {
        std::vector<VertexPtr> neighbor;
        VertexPtr vi = std::make_shared<Vertex>(V[i]);
        neighbors(neighbor, vi);
        for (size_t j = 0; j < neighbor.size(); j++) {
            if (vi->getColor() == neighbor[j]->getColor()) {
                result.emplace_back(vi->getNumber(), neighbor[j]->getNumber());
            }
        }
    }
     */
    // std::cout << "the number of error color is " << result.size();
    // std::cout << std::endl;

    return palette.getNumberOfPigment();
}

int Graph::visitVertexColor(size_t index) const {
    return V[index].getColor();
}

int Graph::visitVertexNumber(size_t index) const {
    return V[index].getNumber();
}

void Graph::parallelColorGraph() {
    std::vector<Vertex> U = V;
}