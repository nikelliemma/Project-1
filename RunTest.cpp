#include <iostream>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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


int main(int argc, char *argv[]){

    char *filename = NULL;
    int option, k = -1, a = -1, L = -1, R = -1, o = -1, S = -1;
    extern char *optarg;
    extern int optopt, optind;

    //GETTING ARGUMENTS FROM COMMAND LINE//

    while((option = getopt(argc, argv, ":f:k:a:l:r:o:s:")) != -1){
        switch (option)
        {
        case 'f':
            if(filename == NULL){
                filename = optarg;
                break;
            }else{
                printf("Argument duplicates. Terminating\n");
                exit(1);
            }
        case 'k':
            if(k == -1){
                if((k = atoi(optarg)) <= 0 ){
                    printf("Malformed input\n");
                    exit(1);
                }
                break;
            }else{
                printf("Argument duplicates. Terminating\n");
                exit(1);
            }
            break;
        case 'a':
            if(a == -1){
                if((a = atoi(optarg)) <= 0 ){
                    printf("Malformed input\n");
                    exit(1);
                }
                break;
            }else{
                printf("Argument duplicates. Terminating\n");
                exit(1);
            }
        case 'l':
            if(L == -1){
                if((L = atoi(optarg)) <= 0 ){
                    printf("Malformed input\n");
                    exit(1);
                }
                break;
            }else{
                printf("Argument duplicates. Terminating\n");
                exit(1);
            }
            break;
        case 'r':
            if(R == -1){
                if((R = atoi(optarg)) <= 0 ){
                    printf("Malformed input\n");
                    exit(1);
                }
                break;
            }else{
                printf("Argument duplicates. Terminating\n");
                exit(1);
            }
            break;
        case 'o':
            if(o == -1){
                o = atoi(optarg);
                if(o < 1 && o > 3){
                    printf("Malformed input\n");
                    exit(1);
                }
                break;
            }else{
                printf("Argument duplicates. Terminating\n");
                exit(1);
            }
            break;
        case 's':
            if(S == -1){
                if(o == 3){
                    if((S = atoi(optarg)) <= 0 ){
                        printf("Malformed input\n");
                        exit(1);
                    }
                    break;
                    
                }else{
                    cout << "Option -s can be used only when option -o is 3" << endl;
                    exit(1);
                }
            }else{
                printf("Argument duplicates. Terminating\n");
                exit(1);
            }
            break;
        case ':':   
            printf("Missing argument. Terminating\n");
            exit(1);
        case '?':
            printf("Unkown argument. Terminating\n");
            exit(1);
        default:
            printf("Unknow error\n");
            exit(1);
        }
    }

    //check for excess arguments
    for(; optind < argc ; optind++){
        printf("Too many arguments. Terminating\n");
        exit(1);
    }


    
    switch(o)
    {
        case 1:
            //run vamana
            break;
        case 2:
            //run filtered
            break;
        case 3:
            if(S == -1){
                printf("Missing argument. Terminating\n");
                exit(1);
            }
            //run stitched
            break;
        default:
            printf("Unknow error. Terminating\n");
            exit(1);
    }
}



/*
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
*/