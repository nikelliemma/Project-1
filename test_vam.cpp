#include <iostream>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>



#include "Vamana.h"

using namespace std;

std::vector<std::vector<int>> readClosestPoints(const std::string& filePath);



int main(int argc, char*argv[]){

    if(argc < 5){
        cout << "Too few arguments!\n";
        return 0;
    }

    int L = std::atoi(argv[1]);
    int R = std::atoi(argv[2]);
    int alpha = std::atoi(argv[3]);
    int k = std::atoi(argv[4]);


    string filename = "dummy-data.bin";

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();

    Vamana v(R,L,alpha);


    auto start = std::chrono::high_resolution_clock::now();

    // RRGraph Vam = v.Filtered_Vamana_Index(d, L, R, alpha);
    // Vam.print_graph();
    // Vam.set_nodes_num(d.get_dataset().size());
    // Vam.write_to_binary_file("vamana_index.bin");

    
    RRGraph Vam1(R);
    Vam1.read_from_binary_file("vamana_index.bin");

    Vam1.print_graph();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    std::cout << "Filtered Vamana index created in : " << duration << " seconds" << std::endl;

    vector<vector<int>> groundtruth = readClosestPoints("Groundtruth.txt");

    map<int, int> filter_map = v.Filtered_Find_Medoid(d.get_dataset(),d.get_filter_set(),1);
    
    string filename_1 = "dummy-queries.bin";

    FilteredDataset q;
    q.set_filepath(filename_1);
    q.read_Query_set();

    vector<Data_Point> queries;
    vector<int> index;
    int counter = 0;

    for(Data_Point point: q.get_dataset()){
        if( point.timestamp == 1 || point.timestamp == 0){
            index.push_back(counter);
        }
        counter++;
    }

    // double recall_sum = 0.0;
    for(int i = 0; i<500;i++){
        LVPair res = v.FilteredGreedySearch(Vam1, filter_map, index[i], k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());
        double recall = v.Get_Recall(groundtruth[index[i]], res.first);
        if(recall >= 75) cout << "recall = " << recall  << endl; 
        //cout << "recall = " << recall  << endl; 
        // recall_sum += recall;
    }        
    // cout << "average recall = " << recall_sum / 500;
    
    return 0;
}








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






















// int main(int argc , char * argv[]){

//     if(argc < 5){
//         cout << "Too few arguments!\n";
//         return 0;
//     }

//     int L = std::atoi(argv[1]);
//     int R = std::atoi(argv[2]);
//     int alpha = std::atoi(argv[3]);
//     int k = std::atoi(argv[4]);

//     Dataset<float> d;
//     d.set_filename("siftsmall_base.fvecs");
//     d.set_type(FLOAT);
//     d.read_dataset();


//     Dataset<int> groundtruth;
//     groundtruth.set_filename("siftsmall_groundtruth.ivecs");
//     groundtruth.set_type(INTEGER);
//     groundtruth.read_dataset();


//     Dataset<float> d1;
//     d1.set_filename("siftsmall_query.fvecs");
//     d1.set_type(FLOAT);
//     d1.read_dataset();


//     auto start = std::chrono::high_resolution_clock::now();

//     Vamana v(R, L, alpha);

//     //v.create_vamana_index("siftsmall_base.fvecs", 50, 20, 1);

//     RRGraph Vam = v.Vamana_Index(d.get_dataset(), L, R, alpha);

//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();



//     //Vam.print_graph();
//     std::cout << "Vamana index created in : " << duration << " seconds" << std::endl;   


//     int medoid = 8736;

//     double avg_recall = 0.0;

//     for(int i=0;i<d1.get_dataset().size(); i++){
//         LVPair res = v.GreedySearch(Vam, medoid, d1.get_vector(i), k, L, d.get_dataset()); 

//         double recall = v.Get_Recall(res.first, groundtruth.get_vector(i));
//         avg_recall += recall;

//         cout << "Query " << i << ": recall = " << recall << "%" << endl;
//     }

//     cout << "average recall : " << (avg_recall / d1.get_dataset().size()) << endl;

//     // for(int i=0;i<100;i++){
//     //     cout << res.first[i] << " ";
//     // }
//     // cout << endl;


//     // for(int i=0;i<100;i++){
//     //     cout << groundtruth.get_data(0,i) << " ";
//     // }
//     // cout << endl;

//     // cout << "VISITED NODES:\n";
//     // for(int node : res.second){
//     //     cout << node << " ";
//     // }
//     // cout << endl; 

//     return 0;
// }

