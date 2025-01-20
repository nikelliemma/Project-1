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

void runStitched(string , int , int ,int ,int , int );
void runVamana( int ,int ,int , int );
void runFiltered( int ,int ,int , int );

std::vector<std::vector<int>> readClosestPoints(const std::string& filePath);

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
            cout << "Running Vamana... " << endl;
            //run vamana
            runVamana(L, R , a, k);
            break;
        case 2:
            cout << "Running Filtered Vamana... " << endl;
            //run filtered
            runFiltered(L, R , a, k);
            break;
        case 3:
            if(S == -1){
                printf("Missing argument. Terminating\n");
                exit(1);
            }
            cout << "Running Stitched Vamana... " << endl;
            runStitched(filename, L, R , a, k, S);
            break;
        default:
            printf("Unknow error. Terminating\n");
            exit(1);
    }

    return 0;
}



















// // Function to compute Euclidean distance between two vectors
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

void runFiltered(int L, int R ,int alpha,int k){
    string filename = "datasets/dummy-data.bin";

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();

    Vamana v(R,L,alpha);

    map<int, int> filter_map = v.Filtered_Find_Medoid(d.get_dataset(),d.get_filter_set(),1);

    string filename_1 = "datasets/dummy-queries.bin";

    FilteredDataset q;
    q.set_filepath(filename_1);
    q.read_Query_set();


    RRGraph Vam(R);
    Vam.read_from_binary_file("filtered_vamana_index.bin");

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

        }
        if(q.get_query_type(i) == 0){

            auto knns = findKNearestNeighborsForQuery(d, q.get_data_point(i), q.get_data_point(i).categorical, k);

            if(knns.size() == 0) continue;

            counter_unfiltered++;

            L = filter_map.size() * 10;
            LVPair res = v.FilteredGreedySearch(Vam, filter_map, i, k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());
            double recall = v.Get_Recall(knns, res.first);

            recall_sum_unfiltered += recall;

        }
    }
    cout << "average recall for filtered queries = " << recall_sum_filtered / counter_filtered << "%" << endl;
    cout << "average recall for unfiltered queries = " << recall_sum_unfiltered / counter_unfiltered << "%" << endl;

    return;
}

void runVamana(int L, int R ,int alpha,int k){

    Dataset<float> d;
    d.set_filename("datasets/siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();


    Dataset<int> groundtruth;
    groundtruth.set_filename("datasets/siftsmall_groundtruth.ivecs");
    groundtruth.set_type(INTEGER);
    groundtruth.read_dataset();


    Dataset<float> d1;
    d1.set_filename("datasets/siftsmall_query.fvecs");
    d1.set_type(FLOAT);
    d1.read_dataset();


    //auto start = std::chrono::high_resolution_clock::now();

    Vamana v(R, L, alpha);

    //v.create_vamana_index("siftsmall_base.fvecs", 50, 20, 1);

    //RRGraph Vam = v.Vamana_Index(d.get_dataset(), L, R, alpha);

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    RRGraph Vam(R);
    Vam.read_from_binary_file("vamana_index.bin");

    //Vam.print_graph();
    //std::cout << "Vamana index created in : " << duration << " seconds" << std::endl;   


    int medoid = 8736;

    double avg_recall = 0.0;

    for(int i=0;i<d1.get_dataset().size(); i++){
        LVPair res = v.GreedySearch(Vam, medoid, d1.get_vector(i), k, L, d.get_dataset()); 

        double recall = v.Get_Recall(res.first, groundtruth.get_vector(i));
        avg_recall += recall;

        cout << "Query " << i << ": recall = " << recall << "%" << endl;
    }

    cout << "average recall : " << (avg_recall / d1.get_dataset().size()) << endl;

    return;
}




void runStitched(string filename, int L, int R ,int alpha,int k, int Rst){

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();

    string filename_1 = "datasets/dummy-queries.bin";

    FilteredDataset q;
    q.set_filepath(filename_1);
    q.read_Query_set();

    Vamana v(R,L,alpha);
    std::vector<std::vector<float>> unfiltered;
    RRGraph Vam(R);
    
    int flag = 0;

    auto start = std::chrono::high_resolution_clock::now();
    
    int size = q.get_dataset().size();
    for(int i = 0; i < size; i++){
        if(q.get_data_point(i).categorical == -1){
            flag++;
        }
    }

    if(flag){   //if there are unfiltered
        size = d.get_dataset().size();
        for(int i = 0; i < size; i++){
            unfiltered.push_back(d.get_data_point(i).data_vector);
        }

        Vam = v.Vamana_Index(unfiltered, L, R, alpha);
        //cout << "Printing Vamana Graph... " << endl;
        //Vam.print_graph();
    }
    
    if(flag < q.get_dataset().size()){ //if there are filtered

        GraphCollection collectionOfGraphs = v.StitchedVamana(d, L, R, Rst, alpha);
        cout << "Printing Stitched Vamana Graphs... " << endl;
        for(int i = 0; i < collectionOfGraphs.size(); i++){
            cout << "Graph for filter " << i << "." << endl;
            collectionOfGraphs[i].print_graph();
        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    std::cout << "Stitched Vamana index created in : " << duration << " seconds" << std::endl;

    /*
    size = d.get_dataset().size();

    for(int i = 0; i < size; i++){

        if(d.get_data_point(i).categorical == -1){

            LVPair res = v.GreedySearch(Vam, find_medoid(unfiltered), d.get_data_point(i).data_vector, k, L, unfiltered);


        }else{

            LVPair res = v.GreedySearch(collectionOfGraphs(d.get_data_point(i).categorical), find_medoid(unfiltered), d.get_data_point(i).data_vector, k, L, unfiltered);
        }
    }
    */
}

