#include <fstream>
#include "../headers/utils.h"

void utils::readCsvData_OneFile(Graph &graph, const std::string &path) {
    std::ifstream stops(path);

    std::string line;
    getline(stops, line);
    while (getline(stops, line)) {
        std::stringstream ss(line);
        std::string originStr, destinyStr, distanceStr;

        getline(ss, originStr, ',');
        getline(ss, destinyStr, ',');
        getline(ss, distanceStr, '\n');

        int origin = std::stoi(originStr);
        int destiny = std::stoi(destinyStr);
        double distance = std::stod(distanceStr);

        graph.addVertex(origin);
        graph.addVertex(destiny);
        graph.addBidirectionalEdge(origin, destiny, distance);
    }
}

void utils::readCsvData_TwoFile(Graph &graph, const std::string &path) {
    std::string e = path + "edges.csv";
    std::string n = path + "nodes.csv";

    std::ifstream nodes(n);
    std::ifstream edges(e);

    std::string line;

    getline(nodes, line);
    getline(edges, line);

    while (getline(nodes, line)) { // NODES
        std::stringstream ss(line);
        std::string idStr, lonStr, latStr;

        getline(ss, idStr, ',');
        getline(ss, lonStr, ',');
        getline(ss, latStr, '\n');

        int id = std::stoi(idStr);
        double lon = std::stod(lonStr);
        double lat = std::stod(latStr);

        graph.addVertex(id, lon, lat);
    }

    while (getline(edges, line)) { // EDGES
        std::stringstream ss(line);
        std::string originStr, destinyStr, distanceStr;

        getline(ss, originStr, ',');
        getline(ss, destinyStr, ',');
        getline(ss, distanceStr, '\n');

        int origin = std::stoi(originStr);
        int destiny = std::stoi(destinyStr);
        double distance = std::stod(distanceStr);

        graph.addBidirectionalEdge(origin, destiny, distance);
    }
}


