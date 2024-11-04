#include <iostream>
#include <chrono>
#include <algorithm>
#include "dataset.h"
#include "Graph.h"

#include "Vamana.h"

using namespace std;



// double get_recall(std::vector<int> vec1, std::vector<int> vec2){

//     //sort both vectors to enable direct comparison
//     // std::sort(vec1.begin(), vec1.end());
//     // std::sort(vec2.begin(), vec2.end());

//     //count matching elements
//     int match_count = 0;
//     for(size_t i = 0; i < vec1.size(); i++){
//         for(int j = 0; j < vec2.size(); j++){
//             if(vec1[i] == vec2[i]){
//                 match_count++;
//             }
//         }
//     }

//     cout << match_count << endl;

//     //calculate percentage of similarity
//     double percentage = (static_cast<double>(match_count) / vec1.size()) * 100;
//     return percentage;
// }




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


    // for(int i = 0; i< d.get_dim(); i++){
    //     cout << d.get_data(0, i) << " ";
    // }
    // cout << endl << endl;

    // cout << "vector dimension : " << d.get_dim() << endl;
    // cout << "total vector num : " << d.get_vectors_num() << endl;

    // RRGraph g(10);
    // g.create_Rregular_graph(d.get_dataset());

    auto start = std::chrono::high_resolution_clock::now();

    Vamana v(50, 50, 1);

    RRGraph Vam = v.Vamana_Index(d.get_dataset(), 100, 50, 1);

    //Vam.print_graph();

    int medoid = 8736;

    LVPair res = v.GreedySearch(Vam, medoid, d1.get_vector(0), 100, 100, d.get_dataset()); 

    double recall = v.Get_Recall(res.first, groundtruth.get_vector(0));

    cout << "recall : " << recall << "%" << endl;

    for(int i=0;i<100;i++){
        cout << res.first[i] << " ";
    }
    cout << endl;


    for(int i=0;i<100;i++){
        cout << groundtruth.get_data(0,i) << " ";
    }
    cout << endl;



    // cout << "VISITED NODES:\n";
    // for(int node : res.second){
    //     cout << node << " ";
    // }
    // cout << endl;


    // LVPair res = v.GreedySearch(g, medoid, 5123, 20, 50, d.get_dataset()); 
    // cout << "NNs of query :" << endl;
    // for(int node : res.first){
    //    cout << node << " ";
    // }
    // cout << endl;
    // cout << endl;

    // cout << "visited nodes :\n\n";
    // for(int node : res.second){
    //    cout << node << " ";
    // }
    // cout << endl;


    // cout << "visited " << res.second.size() << " nodes" << endl;
    // std::cout << std::flush;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Vamana index created in : " << duration << " seconds" << std::endl;




    // v.RobustPruning(g, 5123, res.second, 1, 10, d.get_dataset());

    //g.print_graph();
    

    return 0;
}

