#include <iostream>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "dataset.h"
#include "Graph.h"

#include "Vamana.h"

using namespace std;


// Function to read the file and store closest points into a vector of vectors
std::vector<std::vector<int>> readClosestPoints(const std::string& filePath) {
    // Vector of vectors to store closest points for each query
    std::vector<std::vector<int>> closestPoints;

    // Open the file
    std::ifstream inputFile(filePath);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open file!" << std::endl;
        return closestPoints; // Return empty vector
    }

    // Read the file line by line
    std::string line;
    while (std::getline(inputFile, line)) {
        // Check if the line contains a query
        if (line.find("Query") != std::string::npos) {
            // Find the next line with closest points
            if (std::getline(inputFile, line)) {
                // Parse the points
                std::istringstream ss(line);
                std::vector<int> points;
                std::string point;

                while (std::getline(ss, point, ',')) {
                    try {
                        points.push_back(std::stoi(point));
                    } catch (const std::invalid_argument& e) {
                        //std::cerr << "Error: Invalid number format in file." << std::endl;
                    }
                }

                // Add points to the main vector
                closestPoints.push_back(points);
            }
        }
    }

    // Close the file
    inputFile.close();
    return closestPoints;
}



int main(int argc, char*argv[]){

    if(argc < 5){
        cout << "Too few arguments!\n";
        return 0;
    }

    int L = std::atoi(argv[1]);
    int R = std::atoi(argv[2]);
    int alpha = std::atoi(argv[3]);
    int k = std::atoi(argv[4]);
    int Rst = std::atoi(argv[5]);

    string filename = "dummy-data.bin";

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();

    string filename_1 = "dummy-queries.bin";

    FilteredDataset q;
    q.set_filepath(filename_1);
    q.read_Query_set();

    Vamana v(R,L,alpha);
    int flag = 0;

    auto start = std::chrono::high_resolution_clock::now();
    
    int size = q.get_dataset().size();
    for(int i = 0; i < size; i++){
        if(q.get_data_point(i).categorical == -1){
            flag++;
        }
    }
    

    std::vector<std::vector<float>> unfiltered;
    RRGraph Vam(R);
    /*
    if(flag){   //if there are unfiltered
        size = d.get_dataset().size();
        for(int i = 0; i < size; i++){
            unfiltered.push_back(d.get_data_point(i).data_vector);
        }
        cout << "doing vamana" << endl;
        Vam = v.Vamana_Index(unfiltered, L, R, alpha);
        cout << "out of vamana" << endl;
        Vam.print_graph();
    }
    */
    if(flag < q.get_dataset().size()){ //if there are filtered
        cout << "about to do stitched" << endl;
        GraphCollection collectionOfGraphs = v.StitchedVamana(d, L, R, Rst, alpha);
        cout << "out of stitched" << endl;

        for(int i = 0; i < collectionOfGraphs.size(); i++){
            cout << "im here " << i << endl;
            collectionOfGraphs[i].print_graph();
        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    std::cout << "Stitched Vamana index created in : " << duration << " seconds" << std::endl;

   //vector<vector<int>> groundtruth = readClosestPoints("Groundtruth.txt");


    return 0;
}


