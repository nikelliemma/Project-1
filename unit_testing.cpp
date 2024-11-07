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


void medoid_test(){
    Dataset<float> d;
    d.set_filename("siftsmall_base.fvecs");
    d.set_type(FLOAT);
    d.read_dataset();
    Vamana v(1,1,1);
    TEST_CHECK(v.find_medoid(d.get_dataset()) == 8736);
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



void R_regular_graph_test(){
    TEST_CHECK(check_graph(10));
    TEST_CHECK(check_graph(4));
    TEST_CHECK(check_graph(0));
}


void greedy_search_test(){
    return;
}


TEST_LIST = {

    {"recall test 1", recall_test_1},
    {"recall test 2", recall_test_2},
    {"recall test 3", recall_test_3},
    {"medoid:", medoid_test},
    {"R-regular graph", R_regular_graph_test},
    {"Greedy Search test", greedy_search_test},
    {0}

};

//------------------------------------------------------//