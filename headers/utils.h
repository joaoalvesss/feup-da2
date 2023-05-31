#ifndef PROJETO_2_UTILS_H
#define PROJETO_2_UTILS_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "graph.h"
#include "vertexEdge.h"

using namespace std;

namespace utils{
        void readCsvData_OneFile(Graph &graph, const std::string& path);
        void readCsvData_TwoFile(Graph &graph, const std::string& path);
}
#endif //PROJETO_2_UTILS_H
