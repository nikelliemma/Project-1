#include <iostream>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <queue>
#include <algorithm>

#include "acutest.h"
#include "Vamana.h"


using namespace std;

//---------------- Euclidean Testing ----------------//

void test_euclidean2(void){

    Vamana myVam(1,1,1);
    double result, euc; 

    std::vector<int> v3 = {1,1}, v4 = {1,1}; 
 
    result = std::sqrt(0), euc = myVam.euclidean_distance(v3, v4);
    TEST_CHECK( euc == result);
}

void test_euclidean3(void){

    Vamana myVam(1,1,1);
    double result, euc; 

    std::vector<int> v1 = {6, 2, 16}, v3 = {1,1}; 

    TEST_EXCEPTION(myVam.euclidean_distance(v1, v3), std::exception);
}



//---------------- Recall Testing ----------------//

void test_recall1(void){

    Vamana myVam(1,1,1);
    std::vector<int> v1 = {1, 2, 3, 4, 5}, v2 = {2, 3, 6, 7};

    double received = myVam.Get_Recall(v1, v2), expected = (2 * 100) / 5;
    TEST_CHECK( received == expected);
}   

void test_recall2(void){

    Vamana myVam(1,1,1);
    
    std::vector<int> v1 = {5, 10, 15, 20, 25};
    std::vector<int> v2 = {5, 10, 15, 20, 25, 30, 35, 40};

    double received = myVam.Get_Recall(v1, v2), expected = (5 * 100) / 5;
    TEST_CHECK(received == expected);
}   

void test_recall3(void){

    Vamana myVam(1,1,1);
    
    std::vector<int> v1 = {1, 2, 3, 4, 5};
    std::vector<int> v2 = {6, 7, 8, 9, 10};

    double received = myVam.Get_Recall(v1, v2), expected = (0 * 100) / 5;
    TEST_CHECK(received == expected);
}   


//---------------- Vamana Index Testing ----------------//

double vam_index(int L, int R, float alpha, int k){

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

    Vamana v(R, L, alpha);

    RRGraph Vam = v.Vamana_Index(d.get_dataset(), L, R, alpha);
    
    int medoid = v.find_medoid(d.get_dataset());

    double recall_counter = 0;

    for(int i=0; i<d1.get_dataset().size(); i++){

        LVPair res = v.GreedySearch(Vam, medoid, d1.get_vector(i), k, L, d.get_dataset()); 

        double recall = v.Get_Recall(res.first, groundtruth.get_vector(i));

        recall_counter += recall;
    }

    return recall_counter / k;
}

void test_index(){

    Dataset<float> d;
    d.set_filename("datasets/siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();


    Vamana v(13, 50, 1);

    RRGraph Vam = v.Vamana_Index(d.get_dataset(), 50, 13, 1);

    int counter = 0;

    for(int i = 0; i < Vam.get_nodes_num(); i++){
        if(Vam.get_node(i)->neighbors.size() <= 14) counter++;
    }

    TEST_CHECK(counter == Vam.get_nodes_num());
    
}


//---------------- Medoid Testing ----------------//

void medoid_test(){

    //siftsmall dataset medoid 
    Dataset<float> d;
    d.set_filename("datasets/siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();
    Vamana v(1,1,1);
    // TEST_CHECK(v.find_medoid(d.get_dataset()) == 8736);

    //single vector test
    std::vector<std::vector<float>> dataset1 = {{1.0, 2.0}};
    TEST_CHECK(v.find_medoid(dataset1) == 0);

    //2 identical vectors test
    std::vector<std::vector<float>> dataset2 = {{1.0, 2.0}, {1.0, 2.0}};
    TEST_CHECK(v.find_medoid(dataset2) == 0 || v.find_medoid(dataset2) == 1);

    //3 2-dimensional vectors
    std::vector<std::vector<float>> dataset3 = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 2.0}};
    TEST_CHECK(v.find_medoid(dataset3) == 1);   

    //4 3-dimensional vectors
    std::vector<std::vector<float>> dataset4 = {{1.0, 2.0, 3.0}, {2.0, 3.0, 4.0}, {4.0, 0.0, 1.0}, {3.0, 3.0, 3.0}};
    TEST_CHECK(v.find_medoid(dataset4) == 3);
}


//---------------- Graph Testing ----------------//

bool check_graph(int R){
    Dataset<float> d;
    d.set_filename("datasets/siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();

    RRGraph graph(R);
    graph.create_Rregular_graph(d.get_dataset());

    int counter = 0;
    for(Node * node : graph.get_graph()){

        std::unordered_set<int> neighbors_set;
        for(int neighbor: node->neighbors){
            neighbors_set.emplace(neighbor);
        }

        if(neighbors_set.size() == R) counter++;
    }

    if(counter == graph.get_nodes_num()) return true;
 
    return false;
}

//test to check if the edges of each node in the Vamana graph are <= R
void Vam_Index_Graph_Test(){

    Dataset<float> d;
    d.set_filename("datasets/siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();


    Vamana v(13, 50, 1);

    RRGraph Vam = v.Vamana_Index(d.get_dataset(), 50, 13, 1);

    int counter = 0;
    for(int i = 0 ; i < Vam.get_nodes_num(); i++){
        if(Vam.get_node(i)->neighbors.size() <= 13) counter++;
    }

    TEST_CHECK(counter==Vam.get_nodes_num());

}

void R_regular_graph_test(){
    TEST_CHECK(check_graph(10));
    TEST_CHECK(check_graph(4));
    TEST_CHECK(check_graph(0));
}


void Filtered_Find_Medoid_Test(){

    Vamana v(13, 100, 1);

    string filename = "datasets/dummy-data.bin";

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();

    std::map<int, int> starting_points = v.Filtered_Find_Medoid(d.get_dataset(), d.get_filter_set(), 1);
    
    bool flag = true;

    for(int i = 0; i < starting_points.size(); i++){
        if(starting_points[i] == -1) flag = false;
    }

    TEST_CHECK(starting_points.size() == d.get_filter_set().size());
    TEST_CHECK(flag);

}

void Filtered_Vamana_Test(){

    RRGraph Vam(13);
    Vam.read_from_binary_file("filtered_vamana_index.bin");

    std::vector<Node *> graph = Vam.get_graph();

    TEST_CHECK(graph.size() == 10000);

}

void Filtered_Vamana_Parallel_Test(){

    RRGraph Vam(13);
    Vam.read_from_binary_file("parallel_filtered_vamana_index.bin");

    std::vector<Node *> graph = Vam.get_graph();

    TEST_CHECK(graph.size() == 10000);
}


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






//tests

void Filtered_Queries_Recall_Test(){

    int R = 13;
    int L =150;
    int a = 1;
    int k = 100;
    Vamana v(R,L,a);

    string filename = "datasets/dummy-data.bin";

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();

    std::map<int, int> filter_map = v.Filtered_Find_Medoid(d.get_dataset(), d.get_filter_set(), 1);

    RRGraph Vam(R);
    Vam.read_from_binary_file("filtered_vamana_index.bin");

    string filename_1 = "datasets/dummy-queries.bin";

    FilteredDataset q;
    q.set_filepath(filename_1);
    q.read_Query_set();
 
    int counter_filtered = 0;
    double recall_sum_filtered = 0.0;
    for(int i = 0; i < 500; i++){

        if(q.get_query_type(i) == 1){

            auto knns = findKNearestNeighborsForQuery(d, q.get_data_point(i), q.get_data_point(i).categorical, k);

            if(knns.size() == 0) continue;
            if(q.get_data_point(i).categorical > filter_map.size()) continue;

            counter_filtered++;
            LVPair res = v.FilteredGreedySearch(Vam, filter_map, i, k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());

            double recall = v.Get_Recall(knns, res.first);


            recall_sum_filtered += recall;

            // cout << "query " << i << " recall = " << recall << endl;

        }
    }
    cout << "average recall for 500 random filtered queries (for k = 100) : " << recall_sum_filtered / counter_filtered << "%" << endl;

    TEST_CHECK((recall_sum_filtered / counter_filtered) >= 95);
}

void Unfiltered_Queries_Recall_Tests(){

    int R = 13;
    int L =150;
    int a = 1;
    int k = 20;
    Vamana v(R,L,a);

    string filename = "datasets/dummy-data.bin";

    FilteredDataset d;
    d.set_filepath(filename);
    d.read_Dataset();

    std::map<int, int> filter_map = v.Filtered_Find_Medoid(d.get_dataset(), d.get_filter_set(), 1);

    RRGraph Vam(R);
    Vam.read_from_binary_file("filtered_vamana_index.bin");

    string filename_1 = "datasets/dummy-queries.bin";

    FilteredDataset q;
    q.set_filepath(filename_1);
    q.read_Query_set();
    
    int counter_unfiltered = 0;
    double recall_sum_unfiltered = 0.0;
    for(int i = 0; i < 20; i++){

        if(q.get_query_type(i) == 0){

            auto knns = findKNearestNeighborsForQuery(d, q.get_data_point(i), q.get_data_point(i).categorical, k);

            if(knns.size() == 0) continue;
            //if(q.get_data_point(i).categorical > filter_map.size()) continue;

            counter_unfiltered++;

            L = filter_map.size() * 10;
            LVPair res = v.FilteredGreedySearch(Vam, filter_map, i, k, L, d.get_filter_set(), d.get_dataset(), q.get_dataset());
            double recall = v.Get_Recall(knns, res.first);


            recall_sum_unfiltered += recall;

        }
    }
    cout << "average recall for 20 random unfiltered queries (for k = 20) : " << recall_sum_unfiltered / counter_unfiltered << "%" << endl;

    TEST_CHECK((recall_sum_unfiltered / counter_unfiltered) >= 30);
}

void Stitched_Vamana_Tests(){

}

// void knns_test_1(){
//     // Create a mock dataset
//     std::vector<Data_Point> points = {
//         {{1.0, 2.0}, 0},
//         {{2.0, 3.0}, 0},
//         {{3.0, 4.0}, 1},
//         {{5.0, 6.0}, 1},
//         {{7.0, 8.0}, 2},
//     };

//     //FilteredDataset dataset;

//     // Query point
//     Data_Point query = {{2.5, 3.5}, 1};

//     // Test case 1: k = 2, queryFilter = -1 (no filter)
//     std::vector<int> result = findKNearestNeighborsForQuery(points, query, -1, 2);
//     TEST_CHECK(result.size() == 2);
//     TEST_CHECK(result[0] == 1); // Closest point
//     TEST_CHECK(result[1] == 2); // Second closest point

//     // Test case 2: k = 1, queryFilter = 1
//     result = findKNearestNeighborsForQuery(points, query, 1, 1);
//     TEST_CHECK(result.size() == 1);
//     TEST_CHECK(result[0] == 2); // Only point in category 1
// }


// void test_findKNearestNeighbors_empty_dataset() {
//     // Create an empty dataset
//     std::vector<Data_Point> points;
//     //FilteredDataset dataset(points);

//     // Query point
//     Data_Point query = {{1.0, 1.0}, 0};

//     // Test case: k = 1, queryFilter = -1
//     std::vector<int> result = findKNearestNeighborsForQuery(points, query, -1, 1);
//     TEST_CHECK(result.empty()); // No points to return
// }

// void test_findKNearestNeighbors_edge_cases() {
//     // Create a mock dataset
//     std::vector<Data_Point> points = {
//         {{1.0, 1.0}, 0},
//         {{2.0, 2.0}, 1},
//     };

//     //FilteredDataset dataset(points);

//     // Query point
//     Data_Point query = {{1.5, 1.5}, 0};

//     // Test case 1: k = 1, queryFilter = -1 (no filter)
//     std::vector<int> result = findKNearestNeighborsForQuery(points, query, -1, 1);
//     TEST_CHECK(result.size() == 1);
//     TEST_CHECK(result[0] == 0); // Closest point to query

//     // Test case 2: k = 2, queryFilter = 0 (only category 0)
//     result = findKNearestNeighborsForQuery(points, query, 0, 2);
//     TEST_CHECK(result.size() == 1);
//     TEST_CHECK(result[0] == 0); // Only one point in category 0

//     // Test case 3: k = 5, more neighbors than points
//     result = findKNearestNeighborsForQuery(points, query, -1, 5);
//     TEST_CHECK(result.size() == 2);
//     TEST_CHECK((result[0] == 0 && result[1] == 1) || (result[0] == 1 && result[1] == 0));
// }



// void Vamana_Recall_Test(){

//     int L = 150;
//     int R = 13;
//     int alpha = 1;
//     int k = 100;

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

//     Vamana v(R, L, alpha);


//     //RRGraph Vam = v.Vamana_Index(d.get_dataset(), L, R, alpha);
//     //Vam.write_to_binary_file("vamana_index.bin");

//     RRGraph Vam(13);
//     Vam.read_from_binary_file("vamana_index.bin");

//     int medoid = 8736;

//     double avg_recall = 0.0;

//     for(int i=0;i<500; i++){

//         LVPair res = v.GreedySearch(Vam, medoid, d1.get_vector(i), k, L, d.get_dataset()); 

//         double recall = v.Get_Recall(res.first, groundtruth.get_vector(i));
//         avg_recall += recall;

//         //cout << "Query " << i << ": recall = " << recall << "%" << endl;

//     }

//     cout << "average recall for siftsmall queries file: " << (avg_recall / 500) << endl;

//     TEST_CHECK((avg_recall / 500) >= 97);

//     return;
// }



//------------------------------------------------------//



TEST_LIST = {
    
    // {"Vamana graph test", Vam_Index_Graph_Test},
    {"Starting points test (filtered find medoid parallel)", Filtered_Find_Medoid_Test},
    {"Recall test for filtered queries", Filtered_Queries_Recall_Test},
    {"Recall test for unfiltered queries", Unfiltered_Queries_Recall_Tests},
    {"medoid tests", medoid_test},
    {"R-regular graph", R_regular_graph_test},
    {"Vamana index recall test 1", test_recall1},
    {"Vamana index recall test 2", test_recall2},
    {"Vamana index recall test 3", test_recall3},
    {"Filtered Vamana creation test", Filtered_Vamana_Test},
    {"Filtered Vamana Parallel creation test", Filtered_Vamana_Parallel_Test},
    //{"Recall for siftsmall queries file test", Vamana_Recall_Test},
    {"euclidean2", test_euclidean2},
    {"euclidean3", test_euclidean3},
    {0}

};

//------------------------------------------------------//

