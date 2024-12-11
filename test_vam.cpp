#include <iostream>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>
//#include "dataset.h"
//#include "Graph.h"

#include "Vamana.h"

using namespace std;



// std::vector<std::vector<int>> parseFile(const std::string& file_path) {
//     // Vector of vectors to store closest points for each query
//     std::vector<std::vector<int>> closest_points;

//     // Open the file
//     std::ifstream file(file_path);
//     if (!file.is_open()) {
//         std::cerr << "Error: Unable to open file!" << std::endl;
//         return closest_points; // Return empty if file can't be opened
//     }

//     std::string line;
//     while (std::getline(file, line)) {
//         std::stringstream ss(line);
//         std::string word;

//         // Check if the line contains a query
//         if (line.find("query") != std::string::npos) {
//             std::vector<int> points;

//             // Skip "closest 100 points are:" and move to the next line
//             std::getline(file, line);
//             std::stringstream points_stream(line);

//             while (std::getline(points_stream, word, ',')) {
//                 // Parse each number and add to the vector
//                 try {
//                     points.push_back(std::stoi(word));
//                 } catch (const std::invalid_argument&) {
//                     // Handle any non-integer values gracefully
//                     continue;
//                 }
//             }

//             // Store the closest points in the main vector
//             closest_points.push_back(points);
//         }
//     }

//     file.close();
//     return closest_points;
// }


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


    string filename = "dummy-data.bin";

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();


    // cout << "vector 9999: \n";
    // for(int i = 0; i< d.get_dimension(); i++){
    //     cout << d.get_data_point(9999).data_vector[i] << " ";
    // }
    // cout << endl;
    // cout << "categorical : " << d.get_data_point(9999).categorical << endl;
    //out << d.get_data_point(0).data_vector.size() << endl;

    // for( int i = 0; i<d.get_dataset().size(); i++){
    //    cout << "vector " << i << " categorical : " << d.get_data_point(i).categorical << endl; 
    // }

    // for(int i = 0; i< 10000; i++){
    //     if(d.get_data_point(i).categorical == -1) cout << "LOLOLOLOL\n";
    // }


    Vamana v(R,L,alpha);


    auto start = std::chrono::high_resolution_clock::now();

    // RRGraph Vam(R);

    // Vam.create_Rregular_empty_graph(d.get_dataset());

    // int medoid = v.find_medoid_f(d.get_dataset());

    // LVPair res = v.FilteredGreedySearch(Vam, medoid, 5432, 10, L, d.get_filter_set(), d.get_dataset());

    // for( int node: res.first){
    //     cout << node << " ";
    // }
    // cout << endl;

    // cout << "VISITED NODES\n";

    // for (const auto& element : res.second) {
    //     std::cout << element << " ";
    // }
    // cout << endl;

    RRGraph Vam = v.Filtered_Vamana_Index(d, L, R, alpha);
    Vam.print_graph();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    std::cout << "Filtered Vamana index created in : " << duration << " seconds" << std::endl;

    vector<vector<int>> groundtruth = readClosestPoints("Groundtruth.txt");

    for(int point: groundtruth[0]){
        cout << point << " ";
    }
    cout << endl;
    cout << endl;

    map<int, int> filter_map = v.Filtered_Find_Medoid(d.get_dataset(),d.get_filter_set(),1);

    int counter = 0;
    
    string filename_1 = "dummy-queries.bin";

    FilteredDataset q;
    q.set_filepath(filename_1);
    q.read_Query_set();


    LVPair res = v.FilteredGreedySearch(Vam, filter_map, 0, k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());



    for(int point: res.first){
        cout << point << " ";
    }
    // cout << endl;
    //     double recall = v.Get_Recall(vec, res.first);
    //     cout << "vector " << counter << ":" << "recall = " << recall << endl;

    //     counter++;
    // }
    

    // cout << endl << "TYPE : "<<q.get_data_point(0).timestamp << endl;
    // cout << endl << "CATEGORICAL  : "<<q.get_data_point(0).categorical << endl;



    return 0;
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
