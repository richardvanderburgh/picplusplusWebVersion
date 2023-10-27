//#include <iostream>
//#include "First_EE.hpp"
#include "init.hpp"

// #include <gtest/gtest.h>

//#include <boost/algorithm/algorithm.hpp>
//#include <nlohmann/json.hpp>


int main(int argc, char* argv[]) {


	Init init;
    double L = atoi(argv[1]);
	int N = atoi(argv[2]);
	int nt = atoi(argv[3]);
	double dt = std::stod(argv[4]);
    int ng = atoi(argv[5]);
    int mode = std::stod(argv[6]);
    double V0 = std::stod(argv[7]);
    int numSpecies = std::stod(argv[8]);
    double amplitude = std::stod(argv[9]);
    double VT1 = std::stod(argv[10]);
    int wp1 = std::stod(argv[11]);
    int qm1 = std::stod(argv[12]);

    //if (argc < 2) {
    //    std::cerr << "Please provide the input file path." << std::endl;
    //    return 1;
    //}

    //std::string inputFile = argv[1]; // Get the input file path from the command-line argument

    //std::unordered_map<std::string, int> arguments; // Map to store parsed arguments

    //std::ifstream file(inputFile); // Open the input file
    //if (!file.is_open()) {
    //    std::cerr << "Failed to open the input file." << std::endl;
    //    return 1;
    //}

    //std::string key, valueStr;
    //while (file >> key >> valueStr) {
    //    int value;
    //    try {
    //        value = std::stoi(valueStr); // Convert the value string to an integer
    //    }
    //    catch (const std::exception& e) {
    //        std::cerr << "Error parsing value for argument: " << key << std::endl;
    //        continue;
    //    }

    //    arguments[key] = value; // Store the parsed argument in the map
    //}

    //// Access the parsed arguments
    //int N = arguments["N"];
    //int nt = arguments["nt"];
    //double dt = arguments["dt"];

    //// Print the parsed arguments
    //std::cout << "N: " << N << std::endl;
    //std::cout << "nt: " << nt << std::endl;
    //std::cout << "dt: " << dt << std::endl;

    //file.close(); // Close the input file

	bool success = init.initialize(L, N, nt, dt, ng, mode, V0, numSpecies, amplitude, VT1, wp1, qm1 );

	return 0;
}