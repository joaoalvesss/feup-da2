#include "../headers/menu.h"
#include <iostream>
#include <chrono>

bool Menu::exitApplication;
Graph Menu::graph;

Menu::Menu(){
    graph = *new Graph();
    exitApplication = false;
    utils::readCsvData_OneFile(graph, "../resources/Toy-Graphs/stadiums.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Toy-Graphs/shipping.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Toy-Graphs/tourism.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_25.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_50.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_75.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_100.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_200.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_300.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_400.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_500.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_600.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_700.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_800.csv");
    // utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_900.csv");
    // utils::readCsvData_TwoFile(graph, "../resources/Real-world Graphs/graph1/");
    // utils::readCsvData_TwoFile(graph, "../resources/Real-world Graphs/graph2/");
    // utils::readCsvData_TwoFile(graph, "../resources/Real-world Graphs/graph3/");
}

void Menu::init() {
    std::atexit(finish);

    while (!exitApplication) {
        int option = showMenu();

        switch (option) {
            case 0:
                exitApplication = true;
                break;

            case 1:
            {
                std::vector<int> path;
                std::vector<bool> visited(graph.getNumVertex(), false);
                double min_cost = std::numeric_limits<double>::max();
                std::vector<int> optimal_path;

                path.push_back(0);
                visited[0] = true;

                auto start_time = std::chrono::high_resolution_clock::now();
                graph.tspBT(path, visited, optimal_path, min_cost, 0.0);
                optimal_path.push_back(0);

                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
                std::cout << "\t> Execution Time: " << duration.count() << " seconds" << std::endl;

                std::cout << "\t> Optimal Path: ";
                for (auto it = optimal_path.rbegin(); it != optimal_path.rend(); ++it)
                    std::cout << *it << " ";

                std::cout << std::endl;
                std::cout << "\t> Minimum Cost: " << min_cost << std::endl;
                graph.resetVisits();
                break;
            }


            case 2:
            {

                std::vector<int> path;
                auto start_time = std::chrono::high_resolution_clock::now();
                double dist = graph.triangularApproximation(path);
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
                std::cout << "\t> Execution Time: " << duration.count() << " seconds" << std::endl;
                std::cout << "\t> Path: ";
                for (int i : path) std::cout << i << " ";
                std::cout << std::endl;
                std::cout << "\t> Minimum Cost: " << dist << std::endl;
                graph.resetVisits();

                break;
            }
            case 3:
            {
                auto start_time = std::chrono::high_resolution_clock::now();
                // Step 1: Build the minimum spanning tree (MST)
                Graph mstGraph = graph.prim(0);

                // Step 2: Find the odd degree vertices and connect them
                std::vector<std::pair<int, int>> mpm = mstGraph.findOddDegreeVerticesAndConnect(mstGraph);

                // Step 3: Add the edges from the MPM to the MST
                mstGraph.addMpmEdgesToMst(mpm, mstGraph);

                // Step 4: Find the Eulerian path in the MST
                std::vector<int> eulerianPath;
                mstGraph.findEulerianPath(0, eulerianPath);

                // Print the Eulerian path
                std::cout << "\t> Eulerian Path: ";
                for (int vertex : eulerianPath) std::cout << vertex << " ";
                std::cout << std::endl;

                // Step 5: Build the Hamiltonian path
                std::vector<int> hamiltonianPath;
                Graph::getHamiltonianPath(eulerianPath, hamiltonianPath);

                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
                std::cout << "\t> Execution Time: " << duration.count() << " seconds" << std::endl;

                // Print the Hamiltonian path
                std::cout << "\t> Hamiltonian Path: ";
                for (int vertex : hamiltonianPath) std::cout << vertex << " ";
                std::cout << std::endl;

                // double cost = calculateTotalDistance(hamiltonianPath);

                break;
            }
            default:
                break;
        }
    }
}


int Menu::showMenu() {
    int choice;

    std::cout << "\n\n";
    std::cout << "\t-------------------------- MAIN MENU ---------------------------\n";
    std::cout << "\t[0] Finish execution and quit\n";
    std::cout << "\t[1] Backtracking Algorithm\n";
    std::cout << "\t[2] Triangular Approximation Heuristic\n";
    std::cout << "\t[3] Other Heuristics (Christofides)\n";
    std::cout << "\t----------------------------------------------------------------\n\n";

    std::cout << "\t> Enter your choice: ";
    std::cin >> choice;

    if (!std::cin) { exit(0);}

    std::cout.flush();
    std::cout << "\n";
    std::cout.flush();

    return choice;
}

void Menu::finish() {
    std::cout << "\t> Finishing execution...\n";
}