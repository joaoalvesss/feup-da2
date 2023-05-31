#include <utility>

#include "../headers/vertexEdge.h"

    /*** Vertex - Stations ***/


Vertex::Vertex(int id, double latitude, double longitude){
    this->id = id;
    this->latitude = latitude;
    this->longitude = longitude;
}

bool Vertex::operator<(Vertex & vertex) const {
    return this->dist < vertex.dist;
}


Edge * Vertex::addEdge(Vertex *destiny, Vertex *origin, double distance) {
    auto newEdge = new Edge(origin, destiny, distance);
    adj.push_back(newEdge);
    destiny->incoming.push_back(newEdge);
    return newEdge;
}

Edge * Vertex::removeEdge(int destID) {
    bool removedEdge = false;
    auto it = adj.begin();
    Edge *res;
    bool found = false;
    while (it != adj.end()) {
        Edge *edge = *it;
        Vertex *dest = edge->getDestiny();
        if (dest->getId() == destID) {
            res = new Edge((*it)->getOrigin(), (*it)->getDestiny(), (*it)->getDistance());
            found = true;
            it = adj.erase(it);
            // Also remove the corresponding edge from the incoming list
            auto it2 = dest->incoming.begin();
            while (it2 != dest->incoming.end()) {
                if ((*it2)->getOrigin()->getId() == id) {
                    it2 = dest->incoming.erase(it2);
                }
                else {
                    it2++;
                }
            }
            delete edge;
            removedEdge = true;
        }
        else {
            it++;
        }
    }
    if (!found) return nullptr;
    return res;
}

std::vector<Edge*> Vertex::getAdj(){
    return this->adj;
}

bool Vertex::isVisited() const {
    return this->visited;
}

Edge *Vertex::getPath() const {
    return this->path;
}

std::vector<Edge *> Vertex::getIncoming() const {
    return this->incoming;
}

void Vertex::setVisited(bool visited) {
    this->visited = visited;
}

void Vertex::setPath(Edge *path) {
    this->path = path;
}

int Vertex::getId() const{
    return id;
}

void Vertex::setId(int id){
    this->id = id;
}

void Vertex::setDist(int dist) {
    this->dist = dist;
}

int Vertex::getDist() {
    return dist;
}
    /*** Edge - Network ***/

Edge::Edge(Vertex *origin, Vertex *destiny, double distance){
    this->origin = origin;
    this->destiny = destiny;
    this->distance = distance;
}

Vertex * Edge::getOrigin() const {
    return this->origin;
}

Vertex * Edge::getDestiny() const {
    return this->destiny;
}

double Edge::getDistance() const {
    return this->distance;
}

Edge *Edge::getReverse() const {
    return this->reverse;
}

void Edge::setReverse(Edge *reverse) {
    this->reverse = reverse;
}

double Vertex::getLatitude() const {
    return latitude;
}
double Vertex::getLongitude() const {
    return longitude;
}
