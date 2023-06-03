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

    /* 4.1 */

    /**
     * @brief backtracking solution for TSP - Time complexity: O(n!)
     * @param path used to store the current path being tested
     * @param visited used to store each nodes have already been visited
     * @param optimal_path used to store the optimal path between recursions
     * @param min_cost the current minimum cost for the visited cycle
     * @param current_cost the current cost of the path stored into the path vector
     */
    void tspBT(std::vector<int>& path, std::vector<bool>& visited, std::vector<int>& optimal_path, double &min_cost, double current_cost);

    /* 4.2 */

    /**
     * @brief simple prim algorithm used to generate the MST starting at a certain node - Time Complexity
     * @param s the starting node of the MST
     * @return the MST in form of an already created graph
     */
    Graph prim(int s);

    /**
     * @brief simple Depth-first search algorithm starting at the given node - Time Complexity: O(V + E)
     * @param id the starting node of the DFS
     * @return a vector with the order of visited nodes
     */
    std::vector<int> dfs(int id);

    /**
     * @brief function used to calculate a triangular approximation cost and an approximate path of the given graph - Time Complexity: O((V + E) log V)
     * @param path vector used to store the the path created by the combination of prim + dfs
     * @return the minimum cost or distance of the stored path
     */
    double triangularApproximation(std::vector<int> &path);

    /**
     * @brief function used to calculate the total distance of a given path (mostly the optimal ones) - Time Complexity: O(V)
     * @param path the path that we want to calculate the minimum cost or distance
     * @return total distance or cost of the given path
     */
    double calculateTotalDistance(const std::vector<int> &path) const;

    /**
     * @brief function used to check if nodes have an edge connection them - Time Complexity: O(D), where D is the average degree of the graph
     * @param v1 origin vertex id
     * @param v2 destiny vertex id
     * @return true if they have an edge connecting them, false otherwise
     */
    bool checkConnectedNodes(int v1, int v2) const;


    /* 4.3 */

    /**
     * @brief it calculates the distance between two points on the Earth's surface using the Haversine formula - Time Complexity: O(1)
     * @param lat1 latitude of the first vertex
     * @param lon1 longitude of the first vertex
     * @param lat2 latitude of the second vertex
     * @param lon2 longitude of the second vertex
     * @return the distance between the two geometrical points
     */
    static double haversine(double lat1, double lon1, double lat2, double lon2);

    /**
     * @brief it constructs a subgraph with vertices of odd degree, finds a minimum weight perfect matching, and combines it with an eulerian circuit to form a hamiltonian path - Time Complexity: O((V + E) log V)
     * @return a vector with the hamiltonian path
     */
    std::vector<int> christofides();

    /**
     * @brief it finds an eulerian circuit in the graph using a stack-based iterative algorithm - Time Complexity: O(V + E)
     * @return a vector with the eulerian path
     */
    std::vector<int> findEulerianCircuit();

    /**
     * @brief it constructs a subgraph with vertices of odd degree and finds a minimum weight perfect matching using a greedy algorithm - Time Complexity: function is O(V^2)
     * @return
     */
    Graph minimumWeightPerfectMatching();

protected:
    std::unordered_map<int, Vertex *> vertexSet;
    int findVertexIdx(const int &id) const;
};

#endif //PROJETO_2_GRAPH_H
