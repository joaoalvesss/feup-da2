#include <fstream>
#include "../headers/utils.h"

void utils::readCsvData_OneFile(Graph &graph, const std::string &path) {
    std::ifstream stops(path);

    std::string line;
    if (!stops.is_open()) {
        std::cout << "Failed to open the file." << std::endl;
        return;
    }
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
    std::string n = path + "nodes.csv";
    std::ifstream nodes(n);
    std::string line1;

    if (!nodes.is_open()) {
        std::cout << "Failed to open the nodes file." << std::endl;
        return;
    }
    getline(nodes, line1);
    while (getline(nodes, line1)) { // NODES
        std::stringstream ss(line1);
        std::string idStr, lonStr, latStr;

        getline(ss, idStr, ',');
        getline(ss, lonStr, ',');
        getline(ss, latStr, '\n');

        int id = std::stoi(idStr);
        double lon = std::stod(lonStr);
        double lat = std::stod(latStr);

        graph.addVertex(id, lon, lat);
    }

    std::string e = path + "edges.csv";
    std::ifstream edges(e);
    std::string line2;

    if (!edges.is_open()) {
        std::cout << "Failed to open the edges file." << std::endl;
        return;
    }
    getline(edges, line2);

    while (getline(edges, line2)) { // EDGES
        std::stringstream ss(line2);
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


