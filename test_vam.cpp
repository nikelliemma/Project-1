#include <iostream>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unistd.h>
#include <vector>
#include <unordered_map>
#include <queue>
#include <cmath>




#include "Vamana.h"
#include "filtered_dataset.h"

using namespace std;



// Function to compute Euclidean distance between two vectors
float euclideanDistance(const std::vector<float>& a, const std::vector<float>& b) {
    float dist = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        dist += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return std::sqrt(dist);
}

// Function to find k-nearest neighbors for a single query (returns indices)
std::vector<int> findKNearestNeighborsForQuery(FilteredDataset& dataset, 
                                               const Data_Point& queryPoint, 
                                               int queryFilter, int k) {
    // Preprocess the dataset by grouping indices by `categorical`
    std::unordered_map<int, std::vector<int>> groupedIndices;
    for (int i = 0; i < dataset.get_vectors_num(); ++i) {
        groupedIndices[dataset.get_data_point(i).categorical].push_back(i);
    }

    std::priority_queue<std::pair<float, int>> pq; // Max-heap for distances and indices

    // Determine dataset subset to compare
    if (queryFilter != -1) {
        // Compare only to points with the same filter
        if (groupedIndices.find(queryFilter) != groupedIndices.end()) {
            for (int idx : groupedIndices[queryFilter]) {
                const auto& datasetPoint = dataset.get_data_point(idx);
                float dist = euclideanDistance(queryPoint.data_vector, datasetPoint.data_vector);
                pq.push({dist, idx});
                if (pq.size() > k) {
                    pq.pop();
                }
            }
        }

    } else {
        // Compare to all points (no filter restriction)
        for (const auto& [filter, indices] : groupedIndices) {
            for (int idx : indices) {
                const auto& datasetPoint = dataset.get_data_point(idx);
                float dist = euclideanDistance(queryPoint.data_vector, datasetPoint.data_vector);
                pq.push({dist, idx});
                if (pq.size() > k) {
                    pq.pop();
                }
            }
        }
    }

    // Collect the k-nearest neighbors (indices only)
    std::vector<int> neighborIndices;
    while (!pq.empty()) {
        neighborIndices.push_back(pq.top().second);
        pq.pop();
    }
    std::reverse(neighborIndices.begin(), neighborIndices.end()); // Optional: Sort by ascending distance

    return neighborIndices;
}



// map<int, int> find_starting_points(RRGraph Vam, map<int, int> filter_map, int query, Vamana v, int L, FilteredDataset d, FilteredDataset q){

//     map<int, int> starting_points;

//     for(int i = 0; i < filter_map.size(); i++){
//         LVPair res = v.FilteredGreedySearch(Vam, filter_map, query, 1, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());
//         starting_points[i] = res.first[0];
//     }

//     return starting_points;
// }


// int countSimilarElements(const std::vector<int>& vec1, const std::vector<int>& vec2) {
//     // Use sets to ensure uniqueness and efficient lookup
//     std::unordered_set<int> set1(vec1.begin(), vec1.end());
//     std::unordered_set<int> set2(vec2.begin(), vec2.end());

//     int count = 0;

//     // Iterate through the first set and check if elements exist in the second set
//     for (const int& elem : set1) {
//         if (set2.find(elem) != set2.end()) {
//             count++;  // Increment count for each common element
//         }
//     }

//     return count;
// }


int main(int argc, char *argv[]){

    int option, k = -1, a = -1, L = -1, R = -1, o = -1, S = -1;
    extern char *optarg;
    extern int optopt, optind;

    //GETTING ARGUMENTS FROM COMMAND LINE//

    while((option = getopt(argc, argv, ":k:a:l:r:")) != -1){
        switch (option)
        {
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

    string filename = "datasets/dummy-data.bin";

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();

    Vamana v(R,L,a);


    auto start = std::chrono::high_resolution_clock::now();

    float alpha = 1.0;

    RRGraph Vam = v.Filtered_Vamana_Index_Parallel(d, L, R, alpha);
    Vam.print_graph();
    Vam.set_nodes_num(d.get_dataset().size());
    //Vam.write_to_binary_file("parallel_filtered_vamana_index.bin");

    
    // RRGraph Vam(R);
    // Vam.read_from_binary_file("vamana_index_2.bin");

    // Vam.print_graph();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    std::cout << "Filtered Vamana index created in : " << duration << " seconds" << std::endl;

    map<int, int> filter_map = v.Filtered_Find_Medoid(d.get_dataset(),d.get_filter_set(),1);

    
    
    string filename_1 = "datasets/dummy-queries.bin";

    FilteredDataset q;
    q.set_filepath(filename_1);
    q.read_Query_set();

    vector<Data_Point> queries;
    vector<int> index;
    int counter = 0;


    std::cout << "Filtered Vamana index created in : " << duration << " seconds" << std::endl;
    int counter_filtered = 0;
    int counter_unfiltered = 0;
    double recall_sum_filtered = 0.0;
    double recall_sum_unfiltered = 0.0;
    for(int i = 0; i < q.get_dataset().size();i++){

        if(q.get_query_type(i) == 1){

            auto knns = findKNearestNeighborsForQuery(d, q.get_data_point(i), q.get_data_point(i).categorical, k);

            if(knns.size() == 0) continue;
            if(q.get_data_point(i).categorical > filter_map.size()) continue;

            counter_filtered++;
            LVPair res = v.FilteredGreedySearch(Vam, filter_map, i, k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());

            double recall = v.Get_Recall(knns, res.first);


            recall_sum_filtered += recall;

            cout << "query " << i << " recall = " << recall << endl;

        }
    }
    std::cout << "Filtered Vamana index created in : " << duration << " seconds" << std::endl;
    //     if(q.get_query_type(i) == 0){

    //         auto knns = findKNearestNeighborsForQuery(d, q.get_data_point(i), q.get_data_point(i).categorical, k);

    //         if(knns.size() == 0) continue;
    //         //if(q.get_data_point(i).categorical > filter_map.size()) continue;

    //         counter_unfiltered++;

    //         L = filter_map.size() * 10;
    //         LVPair res = v.FilteredGreedySearch(Vam, filter_map, i, k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());
    //         double recall = v.Get_Recall(knns, res.first);


    //         recall_sum_unfiltered += recall;

    //         //LVPair res = v.FilteredGreedySearch(Vam1, filter_map, i, k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());
    //         //L = filter_map.size() * 10;
    //         //map<int, int> starting_points = find_starting_points(Vam, filter_map, i, v, L, d, q);
    //         //LVPair res = v.FilteredGreedySearch(Vam, starting_points, i, k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());

    //     }
    // }
    // cout << counter_unfiltered << endl;
    // cout << "average recall for filtered queries = " << recall_sum_filtered / counter_filtered << "%" << endl;
    // cout << "average recall for filtered queries = " << recall_sum_unfiltered / counter_unfiltered << "%" << endl;










    //RUN STITCHED VAMANA 

    // start = std::chrono::high_resolution_clock::now();

    // cout << endl << "RUNNING STITCHED VAMANA \n\n";
    // //vector<RRGraph> graphs = v.StitchedVamana(d, 100, 32, 13, 1);
    // GraphCollection graphs = v.StitchedVamanaParallel(d,100, 13, 13, 1);

    // end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    // std::cout << "Stitched Vamana index created in : " << duration << " seconds" << std::endl;

    // vector<vector<float>> dataset2;
    // vector<Data_Point> dataset;
    // for(int i = 0; i < d.get_dataset().size(); i++) {
    //     if(d.get_data_point(i).categorical == 0) {
    //         dataset2.push_back(d.get_data_point(i).data_vector);
    //         dataset.push_back(d.get_data_point(i));
    //     }
    // }

    // RRGraph graph = graphs[0];

    // graph.print_graph();

    // // auto knns = findKNearestNeighborsForQuery(d, q.get_data_point(9980), q.get_data_point(9980).categorical, k);
    // // LVPair res = v.GreedySearch(graph, 1487, q.get_data_point(9980).data_vector, k, L, dataset2);
    // // cout << endl;
    // // for(int point: res.first){
    // //     cout << point << " ";
    // // }
    // // cout << endl;

    // // double recall = v.Get_Recall(knns, res.first);
    // // cout << "query " << 9980 << " recall = " << recall << endl;

    // // int medoid = v.find_medoid_filtered(d.get_dataset(), 0);
    // // int medoid = v.find_medoid(dataset2);
    // // cout << "medoid = " << medoid << endl;

    // cout << "dataset size = " << dataset2.size() << endl;

    // for(int i = 0; i < q.get_dataset().size(); i++){
    //     if(q.get_data_point(i).timestamp != 1) continue;
    //     if(q.get_data_point(i).categorical == 0){
    //         auto knns = findKNearestNeighborsForQuery(d, q.get_data_point(i), q.get_data_point(i).categorical, k);

    //         LVPair res = v.GreedySearch(graph, 1487, q.get_data_point(i).data_vector, k, L, dataset2);

    //         double recall = v.Get_Recall(knns, res.first);
    //         cout << "query " << i << " recall = " << recall << endl;
    //     }
    // }

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
//     d.set_filename("datasets/siftsmall_base.fvecs");
//     d.set_type(FLOAT);
//     d.read_dataset();


//     Dataset<int> groundtruth;
//     groundtruth.set_filename("datasets/siftsmall_groundtruth.ivecs");
//     groundtruth.set_type(INTEGER);
//     groundtruth.read_dataset();


//     Dataset<float> d1;
//     d1.set_filename("datasets/siftsmall_query.fvecs");
//     d1.set_type(FLOAT);
//     d1.read_dataset();


//     auto start = std::chrono::high_resolution_clock::now();

//     Vamana v(R, L, alpha);

//     //v.create_vamana_index("siftsmall_base.fvecs", 50, 20, 1);

//     RRGraph Vam = v.Vamana_Index(d.get_dataset(), L, R, alpha);
//     Vam.write_to_binary_file("vamana_index.bin");

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

