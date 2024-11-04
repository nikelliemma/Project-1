#include <iostream>
#include <chrono>
#include <algorithm>
#include "dataset.h"
#include "Graph.h"

#include "Vamana.h"

using namespace std;


int main(void){


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


    auto start = std::chrono::high_resolution_clock::now();

    Vamana v(50, 50, 1);

    RRGraph Vam = v.Vamana_Index(d.get_dataset(), 100, 50, 1);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Vamana index created in : " << duration << " seconds" << std::endl;   

    //Vam.print_graph();

    int medoid = 8736;

    LVPair res = v.GreedySearch(Vam, medoid, d1.get_vector(0), 50, 50, d.get_dataset()); 

    double recall = v.Get_Recall(res.first, groundtruth.get_vector(0));

    cout << "recall : " << recall << "%" << endl;

    // for(int i=0;i<100;i++){
    //     cout << res.first[i] << " ";
    // }
    // cout << endl;


    // for(int i=0;i<100;i++){
    //     cout << groundtruth.get_data(0,i) << " ";
    // }
    // cout << endl;



    // cout << "VISITED NODES:\n";
    // for(int node : res.second){
    //     cout << node << " ";
    // }
    // cout << endl; 

    return 0;
}

