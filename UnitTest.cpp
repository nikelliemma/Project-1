
#include <iostream>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <stdexcept>

#include "acutest.h"
#include "Vamana.h"


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


//tests for greedy
void test_greedy(void){

}



//tests for pruning
void test_pruning(void){

    RRGraph G(5);
    G.get_node(0)->neighbors = {1, 2, 3, 4};

    std::vector<std::vector<double>> data = {
        {0.0, 0.0}, 
        {1.0, 1.0},
        {2.0, 2.0},
        {3.0, 3.0},
        {4.0, 4.0}
    };

    std::unordered_set<int> V = {1, 2, 3, 4};
    int q = 0, R = 2;
    float a = 1.5;
    Vamana myVam(1,1,1);

    const std::vector<int>& pruned_neighboors = G.get_node(q)->neighbors;

    myVam.RobustPruning(G, q, V, a, R, data);

    TEST_CHECK(((G.get_node(q)->neighbors).size()) <= R);
    
}


TEST_LIST = {
    {"euclidean1", test_euclidean1},
    {"euclidean2", test_euclidean2},
    {"euclidean3", test_euclidean3},
    {"euclidean4", test_euclidean4},
    { "pruning", test_pruning },
    { "greedy", test_greedy },
    {"recall1", test_recall1},  
    {"recall2", test_recall2},
    {"recall3", test_recall3},
    { NULL, NULL }
};