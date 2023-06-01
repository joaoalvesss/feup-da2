#ifndef PROJETO_2_VERTEXEDGE_H
#define PROJETO_2_VERTEXEDGE_H
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <string>

class Edge;

    /*** Vertex ***/

class Vertex {
public:
    Vertex(int id, double latitude = 0.0, double longitude = 0.0);
    bool operator<(Vertex & vertex) const;

    std::unordered_map<int, Edge *> getAdj();
    bool isVisited() const;
    Edge *getPath() const;
    std::unordered_map<int, Edge *> getIncoming() const;
    int getId() const;

    void setId(int id);
    void setVisited(bool visited);
    void setPath(Edge *path);
    void setDist(int dist);
    Edge *addEdge(Vertex *destiny, Vertex *origin, double distance);
    Edge *removeEdge(int destID);
    int getDist();
    int queueIndex = 0;
    double getLatitude() const;
    double getLongitude() const;

protected:
    int id;
    double latitude;
    double longitude;

    // outgoing and coming edges
    std::unordered_map<int, Edge *> adj;
    std::unordered_map<int, Edge *> incoming;

    // auxiliary fields
    bool visited = false;
    Edge *path = nullptr;
    int dist;
};

    /*** Edge ***/

class Edge {
public:
    Edge(Vertex *origin, Vertex *destiny, double distance);

    Vertex * getDestiny() const;
    double getDistance() const;
    Vertex * getOrigin() const;
    Edge *getReverse() const;
    void setReverse(Edge *reverse);

    void getDest();

protected:
    Vertex * destiny;
    double distance;

    // used for bidirectional edges
    Vertex *origin;
    Edge *reverse = nullptr;
};

#endif //PROJETO_2_VERTEXEDGE_H
