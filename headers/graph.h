#ifndef PROJETO_2_GRAPH_H
#define PROJETO_2_GRAPH_H
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <map>
#include <list>
#include <unordered_map>
#include "VertexEdge.h"
#include "MutablePriorityQueue.h"
#include <chrono>
#include <unordered_set>
#include <stack>
#include <set>

class Graph {
public:
    Graph();
    Vertex *findVertex(const int &id) const;
    void addVertex(const int &id, const double& lon = 0.0, const double& lat = 0.0);
    bool removeVertex(const int &id);
    void addEdge(const int &source, const int &dest, double distance) const;
    void addBidirectionalEdge(const int &source, const int &dest, double distance) const;
    int getNumVertex() const;
    std::unordered_map<int, Vertex *> getVertexSet() const;
    void resetVisits();

    /** 4.1 **/
    void tspBT(std::vector<int>& path, std::vector<bool>& visited, std::vector<int>& optimal_path, double &min_cost, double current_cost); // 4.1 DONE
    // std::vector<Vertex *> Graph::prim(Vertex *start);

    /** 4.2 **/
    Graph prim(int s);
    std::vector<int> dfs(int id);
    double triangularApproximation(std::vector<int> &path);
    double calculateTotalDistance(const std::vector<int> &path) const;
    bool check_if_nodes_are_connected(int v1, int v2) const;


    /** 4.3 **/
    static double haversine(double lat1, double lon1, double lat2, double lon2);
    std::vector<int> christofides();
    std::vector<int> findEulerianCircuit();
    Graph minimumWeightPerfectMatching();

protected:
    std::unordered_map<int, Vertex *> vertexSet;
    int findVertexIdx(const int &id) const;
};

#endif //PROJETO_2_GRAPH_H
