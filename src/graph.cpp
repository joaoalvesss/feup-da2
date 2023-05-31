#include "../headers/graph.h"
#include "../headers/utils.h"
#include "../headers/MutablePriorityQueue.h"
#include <cmath>

Graph::Graph(){}

bool Graph::removeVertex(const int &id){
    int idx = findVertexIdx(id);
    if (idx == -1) return false;
    Vertex *vertex = vertexSet[idx];

    for (auto edge : vertex->getAdj()) {
        Vertex *destVertex = edge->getDestiny();
        destVertex->removeEdge(id);
    }
    for (auto edge : vertex->getIncoming()) {
        Vertex *origVertex = edge->getOrigin();
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

void Graph::addEdge(const int &source, const int &dest, double distance) {
    auto v1 = findVertex(source);
    auto v2 = findVertex(dest);
    v1->addEdge(v2, v1, distance);
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
            for (auto edge : vertexSet[last_vertex]->getAdj()) {
                if (edge->getDestiny()->getId() == start_vertex) {
                    double cycle_cost = current_cost + edge->getDistance();
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
        for (auto edge : vertexSet[last_vertex]->getAdj()) {
            if (!visited[edge->getDestiny()->getId()]) {
                path.push_back(edge->getDestiny()->getId());
                visited[edge->getDestiny()->getId()] = true;
                tspBT(path, visited, optimal_path, min_cost, current_cost + edge->getDistance());
                path.pop_back();
                visited[edge->getDestiny()->getId()] = false;
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
        for (auto e : vl->getAdj()) {
            if (!e->getDestiny()->isVisited()) {
                if (e->getDestiny()->getDist() > vl->getDist() + e->getDistance()) {
                    e->getDestiny()->setDist(vl->getDist() + e->getDistance());
                    e->getDestiny()->setPath(e);
                    q.insert(e->getDestiny());
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
            mst.addEdge(i.second->getPath()->getOrigin()->getId(), i.second->getId(), 0);
        }
    }
    return mst;
}


std::vector<int> Graph::dfs(int id) {
    std::vector<int> res;
    Vertex *src = findVertex(id);
    res.push_back(id);
    src->setVisited(true);
    for (auto edg : src->getAdj()) {
        if (!edg->getDestiny()->isVisited()) {
            std::vector<int> a = dfs(edg->getDestiny()->getId());
            for (int i : a) {
                res.push_back(i);
            }
        }
    }
    return res;


    /*std::vector<int> res;
    Vertex *src = findVertex(id);
    res.push_back(id);
    src->setVisited(true);
    std::queue<Vertex *> q;
    q.push(src);
    while (!q.empty()) {
        Vertex *v = q.front();
        q.pop();
        for (auto e : v->getAdj()) {
            if (!e->getDestiny()->isVisited()) {
                q.push(e->getDestiny());
                res.push_back(e->getDestiny()->getId());
                e->getDestiny()->setVisited(true);
            }
        }
    }
    return res;*/
}

double Graph::triangularApproximation(std::vector<int> &path) {
    // Create the MST using Prim's algorithm
    Vertex *parent = findVertex(0);
    Graph mst = prim(0);

    // Perform DFS traversal to obtain the order of visited cities
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

        for (auto edge : v1->getAdj()) {
            if (edge->getDestiny() == v2) {
                totalDistance += edge->getDistance();
                break;
            }
        }
    }

    return totalDistance;
}

bool Graph::check_if_nodes_are_connected(int v1, int v2) const{
    for(const auto& edge : findVertex(v1)->getAdj()){
        if(edge->getDestiny()->getId() == v2)
            return true;
        }
    return false;
}

    /************************************** 4.3 ***************************************/
double Graph::getDistance(int v1, int v2) const {
    for(auto &edge : vertexSet.find(v1)->second->getAdj()) {
        if(edge->getDestiny()->getId() == v2) {
            return edge->getDistance();
        }
    }
    return 0;
}

bool Graph::vertexExists(int vertexID){
    for (const auto& vertex : vertexSet) {
        if (vertex.second->getId() == vertexID) {
            return true;
        }
    }
    return false;
}

std::vector<int> Graph::findOddDegreeVertices(){

   std::vector<int> oddDegreeVertices;

   for (int vertex = 0; vertex < vertexSet.size(); ++vertex) {
       if (vertexSet[vertex]->getAdj().size() % 2 != 0) {
           oddDegreeVertices.push_back(vertex);
       }
   }

   return oddDegreeVertices;
}

void Graph::findEulerianPath(int start_vertex, std::vector<int> &circuit){
    std::stack<int> stack;
    stack.push(start_vertex);

    while (!stack.empty()) {
        int current_vertex = stack.top();

        for(auto edge : findVertex(current_vertex)->getAdj()){
            std::cout <<  "\t" << edge->getOrigin()->getId() << " -> " << edge->getDestiny()->getId() << ": dist -> " << edge->getDistance() << endl;
        }

        if (!vertexSet[current_vertex]->getAdj().empty()) {
            int next_vertex = vertexSet[current_vertex]->getAdj().back()->getDestiny()->getId();
            vertexSet[current_vertex]->getAdj().pop_back();
            stack.push(next_vertex);
        }
        else {
            circuit.push_back(current_vertex);
            stack.pop();
        }
        cout << "\t> CIRCUIT: " << circuit.size() << endl;
        cout << "\t> STACK: " << stack.size() << endl;
    }
}


void Graph::buildMstGraph(Graph &mstGraph, const std::vector<std::pair<int, int>>& mst) const{
    for(auto & i : mst){
        int v1 = i.first;
        int v2 = i.second;

        if(v1 == -1 || v2 == -1) continue;

        if(!mstGraph.vertexExists(v1)){
            mstGraph.addVertex(v1, findVertex(v1)->getLatitude(), findVertex(v1)->getLongitude());

            if(!mstGraph.vertexExists(v2)){
                mstGraph.addVertex(v2, findVertex(v2)->getLatitude(), findVertex(v2)->getLongitude());
                if(!mstGraph.check_if_nodes_are_connected(v1, v2)){
                    double distance = haversine(findVertex(v1)->getLatitude(), findVertex(v1)->getLongitude(), findVertex(v2)->getLatitude(), findVertex(v2)->getLongitude());
                    if(distance == 0.0) {
                        mstGraph.addEdge(v1, v2, getDistance(v1, v2));
                        continue;
                    }
                    mstGraph.addEdge(v1, v2, distance);
                }
                else{
                    mstGraph.addEdge(v1, v2, getDistance(v1, v2));
                }
            }
            else{
                if(!mstGraph.check_if_nodes_are_connected(v1, v2)){
                    double distance = haversine(findVertex(v1)->getLatitude(), findVertex(v1)->getLongitude(), findVertex(v2)->getLatitude(), findVertex(v2)->getLongitude());
                    if(distance == 0.0) {
                        mstGraph.addEdge(v1, v2, getDistance(v1, v2));
                        continue;
                    }
                    mstGraph.addEdge(v1, v2, distance);
                }
                else{
                    mstGraph.addEdge(v1, v2, getDistance(v1, v2));
                }
            }
        }
        else if(!mstGraph.vertexExists(v2)){
            mstGraph.addVertex(v2, findVertex(v2)->getLatitude(), findVertex(v2)->getLongitude());
            if(!mstGraph.check_if_nodes_are_connected(v1, v2)){
                double distance = haversine(findVertex(v1)->getLatitude(), findVertex(v1)->getLongitude(), findVertex(v2)->getLatitude(), findVertex(v2)->getLongitude());
                if(distance == 0.0) {
                    mstGraph.addEdge(v1, v2, getDistance(v1, v2));
                    continue;
                }
                mstGraph.addEdge(v1, v2, distance);
            }
            else{
                mstGraph.addEdge(v1, v2, getDistance(v1, v2));
            }
        }
        else{
            if(!mstGraph.check_if_nodes_are_connected(v1, v2)){
                double distance = haversine(findVertex(v1)->getLatitude(), findVertex(v1)->getLongitude(), findVertex(v2)->getLatitude(), findVertex(v2)->getLongitude());
                if(distance == 0.0) {
                    mstGraph.addEdge(v1, v2, getDistance(v1, v2));
                    continue;
                }
                mstGraph.addEdge(v1, v2, distance);
            }
            else{
                mstGraph.addEdge(v1, v2, getDistance(v1, v2));
            }
        }
    }
}

std::vector<std::pair<int, int>> Graph::findOddDegreeVerticesAndConnect(Graph &mstGraph) const{
    std::vector<int> oddDegreeVertices = mstGraph.findOddDegreeVertices();
    std::vector<std::pair<int, int>> mpm;
    std::vector<bool> visited(oddDegreeVertices.size(), false);

    for(int i = 0; i < oddDegreeVertices.size(); i++){
        double min_distance = std::numeric_limits<double>::max();
        int nearest_neighbor = -1;
        std::pair<int, int> edge; //stores the index of the vertices in oddDegreeVertices

        if(!visited[i]){
            for(int j = 0; j < oddDegreeVertices.size(); j++){
                if(i != j && !visited[j]){
                    double distance = getDistance(oddDegreeVertices[i], oddDegreeVertices[j]);
                    if(distance < min_distance){
                        min_distance = distance;
                        nearest_neighbor = j;
                        edge = std::make_pair(i, nearest_neighbor);
                    }
                }
            }

            if(mstGraph.getDistance(oddDegreeVertices[edge.first], oddDegreeVertices[edge.second]) == 0){
                mpm.emplace_back(oddDegreeVertices[edge.first], oddDegreeVertices[edge.second]);
                visited[edge.first] = true;
                visited[edge.second] = true;
            }

        }
    }

    return mpm;
}

void Graph::addMpmEdgesToMst(const std::vector<std::pair<int, int>>& mpm, Graph &mstGraph) const{
    for(auto & i : mpm){
        double distance = haversine(findVertex(i.first)->getLatitude(), findVertex(i.first)->getLongitude(), findVertex(i.second)->getLatitude(), findVertex(i.second)->getLatitude());

        if(distance == 0.0) mstGraph.addEdge(i.first, i.second, getDistance(i.first, i.second));
        else mstGraph.addEdge(i.first, i.second, distance);
    }
}

void Graph::getHamiltonianPath(const std::vector<int>& eulerian_path, std::vector<int> &hamiltonian_path){
    for (int vertex : eulerian_path) {
        auto it = std::find(hamiltonian_path.begin(), hamiltonian_path.end(), vertex);
        if (it == hamiltonian_path.end()) hamiltonian_path.push_back(vertex);

    }
    hamiltonian_path.push_back(hamiltonian_path.front());
}

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

