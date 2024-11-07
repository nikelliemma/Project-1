#include <iostream>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <stdexcept>

#include "acutest.h"
#include "Vamana.h"




//---------------- Vamana Index Testing ----------------//


double vam_index(int L, int R, float alpha, int k){

    Dataset<float> d;
    d.set_filename("siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();


    Dataset<int> groundtruth;
    groundtruth.set_filename("siftsmall_groundtruth.ivecs");
    groundtruth.set_type(INTEGER);
    groundtruth.read_dataset();


    Dataset<float> d1;
    d1.set_filename("siftsmall_query.fvecs");
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


void recall_test_1(){
    TEST_CHECK(vam_index(100, 13, 1, 100) >= 95);
}

void recall_test_2(){
    TEST_CHECK(vam_index(50, 13, 1, 20) >= 95);
}

void recall_test_3(){
    TEST_CHECK(vam_index(100, 13, 1, 50) >= 95);
}





void test_index(){

    Dataset<float> d;
    d.set_filename("siftsmall_base.fvecs");
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

//test to check if the edges of each node in the Vamana graph are <= R
void Vam_Index_Graph_Test(){

    Dataset<float> d;
    d.set_filename("siftsmall_base.fvecs");
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


void medoid_test(){

    //siftsmall dataset medoid 
    Dataset<float> d;
    d.set_filename("siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();
    Vamana v(1,1,1);
    TEST_CHECK(v.find_medoid(d.get_dataset()) == 8736);

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


bool check_graph(int R){


    Dataset<float> d;
    d.set_filename("siftsmall_base.fvecs");
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

//test to check if every node in the random R-regular graph has R edges
void R_regular_graph_test(){
    TEST_CHECK(check_graph(10));
    TEST_CHECK(check_graph(4));
    TEST_CHECK(check_graph(0));
}

// tests for euclidean
void test_euclidean1(void){

    Vamana myVam(1,1,1);
    double result, euc; 

    std::vector<int> v1 = {6, 2, 16}, v2 = {4, 24, 5};
    
    euc = myVam.euclidean_distance(v1, v2);
    result = std::sqrt(609);

    TEST_CHECK( euc == result);
}

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

void test_euclidean4(void){

    Vamana myVam(1,1,1);
    double result, euc, tolerance = 1e-6; 

    std::vector<float> vf1 = {4.0, 5.67, -14.5}, vf2 = {1.02, -7.5, 9.0};

    result = std::sqrt((4.0 - 1.02) * (4.0 - 1.02) + (5.67 + 7.5) * (5.67 + 7.5) + (-14.5 - 9.0) * (-14.5 - 9.0));
    euc = myVam.euclidean_distance(vf1, vf2);
    TEST_CHECK(std::abs(euc - result) < tolerance);
}


//tests for recall
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


//tests for pruning
void test_pruning(void){

    Dataset<float> d;
    d.set_filename("siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();

    RRGraph graph(50);
    graph.create_Rregular_graph(d.get_dataset());

    Vamana v(13, 50, 1);
    RRGraph Vam = v.Vamana_Index(d.get_dataset(), 50, 13, 1);
    
    //cout <<  endl << Vam.get_graph().size() << endl << graph.get_nodes_num() << endl;
    TEST_CHECK(Vam.get_graph().size() <= graph.get_graph().size());
}



TEST_LIST = {

    {"Vamana index recall test 1", recall_test_1},
    {"Vamana index recall test 2", recall_test_2},
    {"Vamana index recall test 3", recall_test_3},
    {"Vamana graph test", Vam_Index_Graph_Test},
    {"medoid tests", medoid_test},
    {"R-regular graph", R_regular_graph_test},
    {"euclidean1", test_euclidean1},
    {"euclidean2", test_euclidean2},
    {"euclidean3", test_euclidean3},
    {"euclidean4", test_euclidean4},
    { "pruning", test_pruning },
    {"recall1", test_recall1},  
    {"recall2", test_recall2},
    {"recall3", test_recall3},
    {0}

};

//------------------------------------------------------//