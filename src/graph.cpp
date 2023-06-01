#include "../headers/graph.h"
#include "../headers/utils.h"
#include <cmath>

Graph::Graph(){}

bool Graph::removeVertex(const int &id){
    int idx = findVertexIdx(id);
    if (idx == -1) return false;
    Vertex *vertex = vertexSet[idx];

    for (auto it : vertex->getAdj()) {
        Vertex *destVertex = it.second->getDestiny();
        destVertex->removeEdge(id);
    }
    for (auto it : vertex->getIncoming()) {
        Vertex *origVertex = it.second->getOrigin();
        origVertex->removeEdge(id);
    }

    vertexSet.erase(idx);
    delete vertex;

    return true;
}

int Graph::getNumVertex() const {
    return (int) vertexSet.size();
}

std::unordered_map<int, Vertex *> Graph::getVertexSet() const {
    return vertexSet;
}

Vertex * Graph::findVertex(const int &id) const {
    for (auto v : vertexSet)
        if (v.second->getId() == id)
            return v.second;
    return nullptr;
}

int Graph::findVertexIdx(const int &id) const {
    for (const auto pair : vertexSet){
        if (pair.second->getId() == id)
            return pair.first;
    }
    return -1;
}

void Graph::addVertex(const int &id, const double& lon, const double& lat) {
    auto * v = new Vertex(id, lon, lat);
    for(auto const &aux : vertexSet){
        if(aux.first == v->getId())
            return;
    }
    vertexSet.insert({id, v});
}

void Graph::addEdge(Vertex * source, Vertex * dest, double distance) {
    source->addEdge(dest, source, distance);
}

void Graph::addBidirectionalEdge(const int &source, const int &dest, double distance) {
    auto v1 = findVertex(source);
    auto v2 = findVertex(dest);
    auto e1 = v1->addEdge(v2, v1, distance);
    auto e2 = v2->addEdge(v1, v2, distance);
    e1->setReverse(e2);
    e2->setReverse(e1);
}

void Graph::resetVisits() {
    for(auto v: vertexSet){
        (*(v.second)).setVisited(false);
    }
}

/************************************** 4.1 ***************************************/
void Graph::tspBT(std::vector<int>& path, std::vector<bool>& visited, std::vector<int>& optimal_path, double& min_cost, double current_cost) {
    if (path.size() == vertexSet.size()) {
        int start_vertex = path.front();
        int last_vertex = path.back();
        for (auto it : vertexSet[last_vertex]->getAdj()) {
            if (it.second->getDestiny()->getId() == start_vertex) {
                double cycle_cost = current_cost + it.second->getDistance();
                if (cycle_cost < min_cost) {
                    min_cost = cycle_cost;
                    optimal_path = path;
                }
                break;
            }
        }
        return;
    }

    int last_vertex = path.back();
    for (auto it : vertexSet[last_vertex]->getAdj()) {
        if (!visited[it.second->getDestiny()->getId()]) {
            path.push_back(it.second->getDestiny()->getId());
            visited[it.second->getDestiny()->getId()] = true;
            tspBT(path, visited, optimal_path, min_cost, current_cost + it.second->getDistance());
            path.pop_back();
            visited[it.second->getDestiny()->getId()] = false;
        }
    }
}

/************************************** 4.2 ***************************************/

Graph Graph::prim(int s) {
    Vertex * start = findVertex(s);
    MutablePriorityQueue<Vertex> q;
    resetVisits();

    for (auto & i : vertexSet) {
        (i.second)->setDist(INT_MAX);
        (i.second)->setPath(nullptr);
        (i.second)->setVisited(false);
    }
    start->setDist(0);
    start->setVisited(true);
    q.insert(start);
    while (!q.empty()) {
        Vertex *vl = q.extractMin();
        vl->setVisited(true);
        for (auto it : vl->getAdj()) {
            if (!it.second->getDestiny()->isVisited()) {
                if (it.second->getDestiny()->getDist() > vl->getDist() + it.second->getDistance()) {
                    it.second->getDestiny()->setDist(vl->getDist() + it.second->getDistance());
                    it.second->getDestiny()->setPath(it.second);
                    q.insert(it.second->getDestiny());
                }
            }
        }
    }

    Graph mst;
    resetVisits();
    mst.addVertex(start->getId(), 0, 0);
    for (auto & i : vertexSet) {
        mst.addVertex(i.second->getId());
    }
    for (auto & i : vertexSet) {
        if (i.second->getPath() != nullptr) {
            mst.addEdge(i.second->getPath()->getOrigin(), i.second, 0);
        }
    }
    return mst;
}


std::vector<int> Graph::dfs(int id) {
    std::vector<int> res;
    Vertex *src = findVertex(id);
    res.push_back(id);
    src->setVisited(true);
    for (auto it : src->getAdj()) {
        if (!it.second->getDestiny()->isVisited()) {
            std::vector<int> a = dfs(it.second->getDestiny()->getId());
            for (int i : a) {
                res.push_back(i);
            }
        }
    }
    return res;
}

double Graph::triangularApproximation(std::vector<int> &path) {
    Vertex *parent = findVertex(0);
    Graph mst = prim(0);
    mst.resetVisits();
    path = mst.dfs(parent->getId());
    path.push_back(parent->getId());
    double total_distance = calculateTotalDistance(path);

    return total_distance;
}

double Graph::calculateTotalDistance(const std::vector<int> &path) {
    double totalDistance = 0.0;

    for (int i = 0; i < path.size() - 1; i++) {
        Vertex *v1 = findVertex(path[i]);
        Vertex *v2 = findVertex(path[i+1]);

        if(!check_if_nodes_are_connected(path[i], path[i+1])){
            totalDistance += haversine(v1->getLatitude(), v1->getLongitude(), v2->getLatitude(), v2->getLongitude());
            continue;
        }

        for (auto it : v1->getAdj()) {
            if (it.second->getDestiny() == v2) {
                totalDistance += it.second->getDistance();
                break;
            }
        }
    }

    return totalDistance;
}

bool Graph::check_if_nodes_are_connected(int v1, int v2) const{
    for(const auto& it : findVertex(v1)->getAdj()){
        if(it.second->getDestiny()->getId() == v2)
            return true;
    }
    return false;
}

/************************************** 4.3 ***************************************/
double Graph::haversine(double lat1, double lon1, double lat2, double lon2) {
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;

    lat1 = (lat1) * M_PI / 180.0;
    lat2 = (lat2) * M_PI / 180.0;

    double a = pow(sin(dLat / 2), 2) + pow(sin(dLon / 2), 2) * cos(lat1) * cos(lat2);
    double rad = 6371;
    double c = 2 * asin(sqrt(a));
    return rad * c;
}

Graph Graph::minimumWeightPerfectMatching() {
    Graph matching;
    std::unordered_set<int> unmatched;
    std::unordered_set<int> matched;

    for (const auto& vertexPair : vertexSet) {
        int id = vertexPair.first;
        unmatched.insert(id);
    }

    while (!unmatched.empty()) {
        int v = *(unmatched.begin());
        unmatched.erase(v);

        double minDistance = std::numeric_limits<double>::max();
        int closestVertex = -1;

        for (const auto& vertexPair : vertexSet) {
            int id = vertexPair.first;
            if (id != v && unmatched.count(id) > 0) {
                Vertex* u = vertexPair.second;
                double distance = haversine(vertexSet[v]->getLatitude(), vertexSet[v]->getLongitude(), u->getLatitude(), u->getLongitude());
                if (distance < minDistance) {
                    minDistance = distance;
                    closestVertex = id;
                }
            }
        }

        if (closestVertex != -1) {
            unmatched.erase(closestVertex);
            matched.insert(v);
            matched.insert(closestVertex);
            Vertex* v1 = findVertex(v);
            Vertex* v2 = findVertex(closestVertex);
            if (v1 != nullptr && v2 != nullptr) {
                matching.addEdge(v1, v2, minDistance);
                cout << "Entrou" << endl;
                cout << matching.vertexSet.size() << endl;
            }
        }
    }
    return matching;
}

// Helper function to find the Eulerian circuit in a graph
std::vector<int> Graph::findEulerianCircuit() {
    std::vector<int> circuit;

    if (vertexSet.empty())
        return circuit;

    std::unordered_map<int, std::vector<int>> adjacencyList;
    for (const auto& vertexPair : vertexSet) {
        int id = vertexPair.first;
        const std::unordered_map<int, Edge*>& adjEdges = vertexPair.second->getAdj();
        std::vector<int> neighbors;
        for (const auto& edge : adjEdges)
            neighbors.push_back(edge.second->getDestiny()->getId());
        adjacencyList[id] = neighbors;
    }

    int startVertex = vertexSet.begin()->first;
    std::stack<int> stack;
    std::vector<int> circuitTemp;
    stack.push(startVertex);

    while (!stack.empty()) {
        int v = stack.top();

        if (!adjacencyList[v].empty()) {
            stack.push(adjacencyList[v].back());
            adjacencyList[v].pop_back();
        } else {
            circuitTemp.push_back(v);
            stack.pop();
        }
    }

    for (auto it = circuitTemp.rbegin(); it != circuitTemp.rend(); ++it)
        circuit.push_back(*it);

    return circuit;
}

// Christofides algorithm to find an approximate Hamiltonian path
std::vector<int> Graph::christofides() {
    std::vector<int> path;

    // Step 1: Create the minimum-weight perfect matching
    Graph matching = minimumWeightPerfectMatching();

    // Step 2: Create a subgraph of odd-degree vertices from the matching
    Graph subgraph;
    std::unordered_set<int> oddDegreeVertices;

    for (const auto& vertexPair : matching.getVertexSet()) {
        int id = vertexPair.first;
        Vertex* vertex = vertexPair.second;
        if (vertex->getAdj().size() % 2 != 0) {
            oddDegreeVertices.insert(id);
            subgraph.addVertex(id, vertex->getLongitude(), vertex->getLongitude());
        }
    }

    for (const auto& vertexPair : matching.getVertexSet()) {
        Vertex* vertex = vertexPair.second;
        const std::unordered_map<int, Edge*>& adjEdges = vertex->getAdj();
        for (const auto& edge : adjEdges) {
            int destId = edge.second->getDestiny()->getId();
            if (oddDegreeVertices.count(destId) > 0) {
                subgraph.addEdge(findVertex(vertexPair.first), edge.second->getDestiny(), edge.second->getDistance());
            }
        }
    }

    // Step 3: Find the Eulerian circuit in the subgraph
    std::vector<int> circuit = subgraph.findEulerianCircuit();

    // Step 4: Create the Hamiltonian path by removing repeated vertices
    std::unordered_set<int> visited;
    for (int vertex : circuit) {
        if (visited.count(vertex) == 0) {
            path.push_back(vertex);
            visited.insert(vertex);
        }
    }

    return path;
}
