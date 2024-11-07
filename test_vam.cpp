#include <iostream>
#include <chrono>
#include <algorithm>
//#include "dataset.h"
//#include "Graph.h"

#include "Vamana.h"

using namespace std;


int main(int argc , char * argv[]){

    if(argc < 5){
        cout << "Too few arguments!\n";
        return 0;
    }

    int L = std::atoi(argv[1]);
    int R = std::atoi(argv[2]);
    int alpha = std::atoi(argv[3]);
    int k = std::atoi(argv[4]);

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

    Vamana v(R, L, alpha);

    //v.create_vamana_index("siftsmall_base.fvecs", 50, 20, 1);

    RRGraph Vam = v.Vamana_Index(d.get_dataset(), L, R, alpha);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();



    Vam.print_graph();
    std::cout << "Vamana index created in : " << duration << " seconds" << std::endl;   


    int medoid = 8736;

    double avg_recall = 0.0;

    for(int i=0;i<d1.get_dataset().size(); i++){
        LVPair res = v.GreedySearch(Vam, medoid, d1.get_vector(i), k, L, d.get_dataset()); 

        double recall = v.Get_Recall(res.first, groundtruth.get_vector(i));
        avg_recall += recall;

        cout << "Query " << i << ": recall = " << recall << "%" << endl;
    }

    cout << "average recall : " << (avg_recall / d1.get_dataset().size()) << endl;

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

