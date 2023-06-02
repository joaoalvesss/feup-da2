#include "../headers/menu.h"
#include <iostream>
#include <chrono>

bool Menu::exitApplication;
Graph Menu::graph;

Menu::Menu(){
    graph = *new Graph();
    exitApplication = false;
    utils::readCsvData_OneFile(graph, "../resources/Toy-Graphs/stadiums.csv");
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
                Graph helper = graph.prim(0);

                std::vector<int> hamiltonianPath = helper.christofides();
                hamiltonianPath.push_back(0);

                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
                std::cout << "\t> Execution Time: " << duration.count() << " seconds" << std::endl;

                std::cout << "\t> Hamiltonian Path: ";
                for (int vertex : hamiltonianPath) std::cout << vertex << " ";
                std::cout << std::endl;

                double cost = graph.calculateTotalDistance(hamiltonianPath);
                std::cout << "\t> Minimum Cost: " << cost << endl;

                break;
            }
            case 4: {
                cout << "\t'''''''''''''''''''''''''''''\n";
                cout << "\t> [0] Go Back\n";
                cout << "\t> [1] Stadiums\n";
                cout << "\t> [2] Shipping\n";
                cout << "\t> [3] Tourism\n";
                cout << "\t> [4] Edges_25\n";
                cout << "\t> [5] Edges_50\n";
                cout << "\t> [6] Edges_75\n";
                cout << "\t> [7] Edges_100\n";
                cout << "\t> [8] Edges_200\n";
                cout << "\t> [9] Edges_300\n";
                cout << "\t> [10] Edges_400\n";
                cout << "\t> [11] Edges_500\n";
                cout << "\t> [12] Edges_600\n";
                cout << "\t> [13] Edges_700\n";
                cout << "\t> [14] Edges_800\n";
                cout << "\t> [15] Edges_900\n";
                cout << "\t> [16] Real-world -> Graph1\n";
                cout << "\t> [17] Real-world -> Graph2\n";
                cout << "\t> [18] Real-world -> Graph3\n";
                cout << "\t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n\n";
                cout << "\t> Enter the number of file: ";
                string path;
                std::cin >> path;
                graph = Graph();

                if(path == "0"){
                    showMenu();
                }
                else if(path == "1"){
                    utils::readCsvData_OneFile(graph, "../resources/Toy-Graphs/stadiums.csv");
                }
                else if(path == "2"){
                    utils::readCsvData_OneFile(graph, "../resources/Toy-Graphs/shipping.csv");
                }
                else if(path == "3"){
                    utils::readCsvData_OneFile(graph, "../resources/Toy-Graphs/tourism.csv");
                }
                else if(path == "4"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_25.csv");
                }
                else if(path == "5"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_50.csv");
                }
                else if(path == "6"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_75.csv");
                }
                else if(path == "7"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_100.csv");
                }
                else if(path == "8"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_200.csv");
                }
                else if(path == "9"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_300.csv");
                }
                else if(path == "10"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_400.csv");
                }
                else if(path == "11"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_500.csv");
                }
                else if(path == "12"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_600.csv");
                }
                else if(path == "13"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_700.csv");
                }
                else if(path == "14"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_800.csv");
                }
                else if(path == "15"){
                    utils::readCsvData_OneFile(graph, "../resources/Extra_Fully_Connected_Graphs/edges_900.csv");
                }
                else if(path == "16"){
                    utils::readCsvData_TwoFile(graph, "../resources/Real-world Graphs/graph1/");
                }
                else if(path == "17"){
                    utils::readCsvData_TwoFile(graph, "../resources/Real-world Graphs/graph2/");
                }
                else if(path == "18"){
                    utils::readCsvData_TwoFile(graph, "../resources/Real-world Graphs/graph3/");
                }


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
    std::cout << "\t[3] Other Heuristic\n";
    std::cout << "\t[4] Change selected graph\n";
    std::cout << "\t----------------------------------------------------------------\n\n";

    std::cout << "\t> Enter your choice: ";
    std::cin >> choice;

    if (choice == '0') { exit(0);}

    std::cout.flush();
    std::cout << "\n";
    std::cout.flush();

    return choice;
}

void Menu::finish() {
    std::cout << "\t> Finishing execution...\n";
}