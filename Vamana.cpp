#include <iostream>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <random>
#include <functional>
#include <queue>
#include <map>
#include <set>

#include "Vamana.h"
#include "Graph.h"
#include "dataset.h"




//constructor
Vamana::Vamana(int R): vamana_index(R){
    this->filename = "";
    this->alpha = 0;
    this->L = 0;
    this->R = R;
    return;
}



//destructor
Vamana::~Vamana(){

}


int Vamana::get_L(){
    return this->L;
}

int Vamana::get_R(){
    return this->R;
}

int Vamana::get_alpha(){
    return this->alpha;
}

void Vamana::set_L(int L){
    this->L = L;
}

void Vamana::set_R(int R){
    this->R = R;
}

void Vamana::set_alpha(int alpha){
    this->alpha = alpha;
}

//templated function to compute Euclidean distance between two vectors
//this function stops earlier if the distance has reached a greater value than the current minimum distance
// template <typename type>
// double euclidean_distance(const std::vector<type>& vec1, const std::vector<type>& vec2, double min_distance){

//     double dist = 0.0;
//     for(int i = 0; i < vec1.size(); ++i){

//         double diff = vec1[i] - vec2[i];
//         dist += diff * diff;

//         if(dist > min_distance) break;

//     }

//     return std::sqrt(dist);
//     //return dist;
// }

//templated function to compute Euclidean distance between two vectors
template <typename type>
double euclidean_distance(const std::vector<type>& vec1, const std::vector<type>& vec2){

    double dist = 0.0;
    for(int i = 0; i < vec1.size(); ++i){

        double diff = vec1[i] - vec2[i];
        dist += diff * diff;

    }

    return std::sqrt(dist);
    //return dist;
}



//templated function to find the medoid of the given dataset of vectors
template <typename type>
int Vamana::find_medoid(std::vector<std::vector<type> > dataset){

    int n = dataset.size();
    int d = dataset[0].size();
    int medoid_index = -1;
    double min_total_distance = std::numeric_limits<double>::max();
    
    for(int i = 0; i < n; ++i){

        double total_distance = 0.0;
        
        for(int j = 0; j < n; ++j){
            if(i != j){
                total_distance += euclidean_distance(dataset[i], dataset[j]);
            }
        }
        
        if(total_distance < min_total_distance){
            min_total_distance = total_distance;
            medoid_index = i;
        }
    }
    
    return medoid_index;    
}


// bool has_unvisited_elements(
//     const std::priority_queue<std::pair<double, int>, 
//                               std::vector<std::pair<double, int>>, 
//                               std::greater<std::pair<double, int>>>& min_heap, 
//     const std::unordered_set<int>& visited) 
// {
//     // Use a temporary min-heap to inspect elements without modifying the original
//     auto temp_heap = min_heap;
    
//     // Limit the number of elements we check to avoid excessive copying
//     const int max_checks = 10; // arbitrary limit, can be tuned
//     int checks = 0;

//     // Check if there's an unvisited element in the temporary heap
//     while (!temp_heap.empty() && checks < max_checks) {
//         int node_id = temp_heap.top().second;
//         temp_heap.pop();
        
//         if (visited.find(node_id) == visited.end()) {
//             return true; // Found an unvisited element
//         }
//         checks++;
//     }

//     return false; // No unvisited elements found within the first few entries
// }

bool has_unvisited_elements(std::priority_queue<std::pair<double, int>, 
                        std::vector<std::pair<double, int>>, 
                        std::greater<std::pair<double, int>>>& min_heap, 
                        const std::unordered_set<int>& visited){
    
    // Create a copy of the min_heap to avoid modifying the original heap
    auto min_heap_copy = min_heap;

    // Iterate over the elements in the min_heap
    while(!min_heap_copy.empty()){

        //get the top element (minimum element in the min-heap)
        int node_id = min_heap_copy.top().second; 
        min_heap_copy.pop(); 

        //check if this element is not in the visited set
        if(visited.find(node_id) == visited.end()){
            return true; // Found an unvisited element
        }
    }

    return false; // All elements in the min_heap are visited
}


template <typename type>
LVPair Vamana::GreedySearch(RRGraph graph, int starting_node, int query, int k, int L, std::vector<std::vector<type> > dataset){

    //result vector contains the NNs of the query
    //visited set contains the visited nodes 
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
    std::unordered_set<int> visited;
    std::vector<int> result;
    std::unordered_set<int> min_heap_nodes;

    min_heap.emplace(euclidean_distance(dataset[starting_node], dataset[query]), starting_node);

    //start the greedy search loop
    while(has_unvisited_elements(min_heap, visited)){

        //get the minimum distance of the nodes in L set that have node been visited
        auto min_heap_copy = min_heap;
        int temp_node;
        while(!min_heap_copy.empty()){

            temp_node = min_heap_copy.top().second;

            if(visited.find(temp_node) != visited.end()){
                min_heap_copy.pop();
            }
            else{ break; }
        }

        int min_idx = temp_node;

        //insert the neighbors of the min distance node to the query into the result (L) set
        for(int node: graph.get_node(min_idx)->neighbors){ 
            if(min_heap_nodes.find(node) == min_heap_nodes.end()){
                double distance = euclidean_distance(dataset[node], dataset[query]);
                min_heap.emplace(distance, node);
                min_heap_nodes.emplace(node);
            }
        }

        //update the visited set of visited nodes
        visited.insert(min_idx);

        //retain the L closest points to the query
        if(min_heap.size() > L){
            auto min_heap_copy = min_heap;
            min_heap = {};
            for(int i = 0; i < L; i++){
                min_heap.emplace(min_heap_copy.top());
                min_heap_copy.pop();
            }    
        }

    }

    //get the k nearest neighbors
    int i = 0;
    while(i < k){
        result.push_back(min_heap.top().second);
        min_heap.pop();
        i++;
    }
    
    //return the NNs and the visited nodes set
    return {result, visited};
}


template <typename type>
LVPair Vamana::GreedySearch(RRGraph graph, int starting_node, std::vector<type> q_vec, int k, int L, std::vector<std::vector<type> > dataset){

    //result vector contains the NNs of the query
    //visited set contains the visited nodes 
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
    std::unordered_set<int> visited;
    std::vector<int> result;
    std::unordered_set<int> min_heap_nodes;

    min_heap.emplace(euclidean_distance(dataset[starting_node], q_vec), starting_node);

    //start the greedy search loop
    while(has_unvisited_elements(min_heap, visited)){

        //get the minimum distance of the nodes in L set that have node been visited
        auto min_heap_copy = min_heap;
        int temp_node;
        while(!min_heap_copy.empty()){

            temp_node = min_heap_copy.top().second;
            if(visited.find(temp_node) != visited.end()){
                min_heap_copy.pop();
            }
            else{ break; }
        }

        int min_idx = temp_node;

        //insert the neighbors of the min distance node to the query into the result (L) set
        for(int node: graph.get_node(min_idx)->neighbors){
            if(min_heap_nodes.find(node) == min_heap_nodes.end()){    
                double distance = euclidean_distance(dataset[node], q_vec);
                min_heap.emplace(distance, node);
                min_heap_nodes.emplace(node);
            }
        }

        //update the visited set of visited nodes
        visited.insert(min_idx);

        //retain the L closest points to the query
        if(min_heap.size() > L){
            auto min_heap_copy = min_heap;
            min_heap = {};
            for(int i = 0; i < L; i++){
                min_heap.emplace(min_heap_copy.top());
                min_heap_copy.pop();
            }    
        }

    }

    //get the k nearest neighbors
    int i = 0;
    while(i < k){
        result.push_back(min_heap.top().second);
        min_heap.pop();
        i++;
    }
    
    //return the NNs and the visited nodes set
    return {result, visited};
}

std::vector<int> get_random_permutation(int n){

    std::vector<int> perm;

    //fill the vector with elements from 0 to n-1
    for(int i = 0; i<n; i++){
        perm.push_back(i);
    }

    std::random_device rd;  //obtain a random seed
    std::mt19937 g(rd());   //seed the generator

    //shuffle the vector
    std::shuffle(perm.begin(), perm.end(), g);

    return perm;

}

template <typename type>
void Vamana::RobustPruning(RRGraph graph, int query, std::unordered_set<int> V_set, float alpha, int R, std::vector<std::vector<type> > dataset){

    std::vector<int> N_out = graph.get_node(query)->neighbors;
    V_set.insert(N_out.begin(), N_out.end());
    V_set.erase(query);

    N_out.clear();

    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;

    while(!V_set.empty()){

        min_heap = {};

        for(int node : V_set){
            double dis = euclidean_distance(dataset[query], dataset[node]);
            min_heap.emplace(dis, node);
        }

        int min_idx = min_heap.top().second;

        if(std::find(N_out.begin(), N_out.end(), min_idx) == N_out.end()){
            N_out.push_back(min_idx);
        }

        if(N_out.size() == R) break;

        std::unordered_set<int> temp_set = V_set;
        for(int node: temp_set){

            double dis_1 = euclidean_distance(dataset[min_idx], dataset[node]);
            double dis_2 = euclidean_distance(dataset[query], dataset[node]);

            if(alpha * dis_1 <= dis_2){
                V_set.erase(node);
            }
        }

    }

    //update the real graph neighbors with the new N_out neighbors   
    graph.get_node(query)->neighbors = std::move(N_out);

    return;
}

template <typename type>
RRGraph Vamana::Vamana_Index(std::vector<std::vector<type> > dataset, int L, int R, float a){

    //create a R-regular graph based on the geometric properties of the dataset
    RRGraph graph(R);
    graph.create_Rregular_graph(dataset);
    //graph.print_graph();

    //get the size of the dataset()
    int N = dataset.size();

    //get the medoid of the dataset 
    // int medoid_index = find_medoid(dataset);
    int medoid_index = 8736;

    //get the random permutation as a starting 
    std::vector<int> perm = std::move(get_random_permutation(N)); 

    //start the Vamana loop 
    for(int i = 0; i < N; i++){

        cout << i << endl;

        //run greedy search 
        LVPair greedy_result = GreedySearch(graph, medoid_index, perm[i], 1, L, dataset);

        //run robust pruning 
        RobustPruning(graph, perm[i], greedy_result.second, a, R, dataset);

        //get the neighbors of perm[i]
        std::vector<int> perm_i_neighbors = graph.get_node(perm[i])->neighbors;
        int n = perm_i_neighbors.size();
        std::vector<int> temp_vec;

        //iterate through all the perm[i] neighbors 
        for(int j: perm_i_neighbors){

            temp_vec = graph.get_node(j)->neighbors;

            if(std::find(temp_vec.begin(), temp_vec.end(), perm[i]) == temp_vec.end()){
                temp_vec.push_back(perm[i]);
            }

            if(temp_vec.size() > R){
                std::unordered_set<int> temp_set(temp_vec.begin(), temp_vec.end());
                RobustPruning(graph, j, temp_set, a, R, dataset);
            }
            else{
                graph.get_node(j)->neighbors.push_back(perm[i]);
            }

        }

    }

    return graph;

}

// std::string extract_format(std::string &file){

//     std::string format;
//     int len = file.length();

//     int pos = file.find_last_of("/");
//     std::string file_string = file.substr(pos + 1);

//     int index = 0;
//     while(file[index] != '.'){ index++; }

//     if(index == len - 1) return "";

//     format = file.substr(index + 1);

//     return format;

// }



// void Vamana::create_vamana_index(std::string filepath){

//     std::string file_format = extract_format(filepath);

//     if(file_format == "fvecs"){
//         Dataset<float> d;
//         d.set_filename(filepath);
//         d.read_dataset();
        
//     }
//     else if(file_format == "ivecs"){
//         Dataset<int> d;
//         d.set_filename(filepath);
//         d.read_dataset();

//     }
//     else if(file_format == "bvecs"){
//         Dataset<unsigned char> d;
//         d.set_filename(filepath);
//         d.read_dataset();

//     }

// }




//explicit instantiation for `float` type
template LVPair Vamana::GreedySearch<float>(RRGraph graph, int starting_node, int query, int k, int L, std::vector<std::vector<float> > dataset);
template LVPair Vamana::GreedySearch<float>(RRGraph graph, int starting_node, std::vector<float> q_vec, int k, int L, std::vector<std::vector<float> > dataset);
template int Vamana::find_medoid<float>(std::vector<std::vector<float>> dataset);
template RRGraph Vamana::Vamana_Index<float>(std::vector<std::vector<float>> dataset, int L, int R, float a);

//explicit instantiation for `int` type
template LVPair Vamana::GreedySearch<int>(RRGraph graph, int starting_node, int query, int k, int L, std::vector<std::vector<int> > dataset);
template LVPair Vamana::GreedySearch<int>(RRGraph graph, int starting_node, std::vector<int> q_vec, int k, int L, std::vector<std::vector<int> > dataset);
template int Vamana::find_medoid<int>(std::vector<std::vector<int> > dataset);
template RRGraph Vamana::Vamana_Index<int>(std::vector<std::vector<int> > dataset, int L, int R, float a);

//explicit instantiation for `unsigned char` type
template LVPair Vamana::GreedySearch<unsigned char>(RRGraph graph, int starting_node, int query, int k, int L, std::vector<std::vector<unsigned char> > dataset);
template LVPair Vamana::GreedySearch<unsigned char>(RRGraph graph, int starting_node, std::vector<unsigned char> q_vec, int k, int L, std::vector<std::vector<unsigned char> > dataset);
template int Vamana::find_medoid<unsigned char>(std::vector<std::vector<unsigned char> > dataset);
template RRGraph Vamana::Vamana_Index<unsigned char>(std::vector<std::vector<unsigned char> > dataset, int L, int R, float a);










// template <typename type>
// void Vamana::RobustPruning(RRGraph graph, int query, std::unordered_set<int> V_set, float alpha, int R, std::vector<std::vector<type> > dataset){

//     std::vector<int> N_out = graph.get_node(query)->neighbors;
//     V_set.insert(N_out.begin(), N_out.end());
//     V_set.erase(query);

//     N_out.clear();

//     std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;

//     while(!V_set.empty()){

//         min_heap = {};

//         for(int node : V_set){
//             double dis = euclidean_distance(dataset[query], dataset[node]);
//             min_heap.emplace(dis, node);
//         }

//         int min_idx = min_heap.top().second;

//         if(std::find(N_out.begin(), N_out.end(), min_idx) == N_out.end()){
//             N_out.push_back(min_idx);
//         }

//         if(N_out.size() == R) break;

//         std::unordered_set<int> temp_set = V_set;
//         for(int node: temp_set){

//             double dis_1 = euclidean_distance(dataset[min_idx], dataset[node]);
//             double dis_2 = euclidean_distance(dataset[query], dataset[node]);

//             if(alpha * dis_1 <= dis_2){
//                 V_set.erase(node);
//             }
//         }

//     }

//     //update the real graph neighbors with the new N_out neighbors 
//     std::unordered_set<int> unique_neighbors(N_out.begin(), N_out.end());
//     N_out.assign(unique_neighbors.begin(), unique_neighbors.end());     

//     graph.get_node(query)->neighbors = N_out;

//     return;
// }



// auto find_in_vector(std::vector<std::pair<int, double>>& distances, int node){
//     auto it = std::find_if(distances.begin(), distances.end(), [node](const std::pair<int, double>& p) {
//         return p.first == node; // Check if the node matches
//     });
//     return it;
// }

// template <typename type>
// void Vamana::RobustPruning(RRGraph graph, int query, std::unordered_set<int> V_set, float alpha, int R, std::vector<std::vector<type> > dataset){

//     std::vector<int> N_out = graph.get_node(query)->neighbors;
//     V_set.insert(N_out.begin(), N_out.end());
//     V_set.erase(query);

//     std::vector<std::pair<int, double>> visited_nodes;

//     for(int node : V_set){
//         double dis = euclidean_distance(dataset[query], dataset[node]);
//         visited_nodes.emplace_back(node, dis);
//     }

//     std::sort(visited_nodes.begin(), visited_nodes.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b){
//         return a.second < b.second; 
//     });

//     N_out.clear();

//     while(!visited_nodes.empty()){

//         int min_idx = visited_nodes[0].first;

//         if(std::find(N_out.begin(), N_out.end(), min_idx) == N_out.end()){
//             N_out.push_back(min_idx);
//         }

//         if(N_out.size() == R) break;

//         std::unordered_set<int> temp_set = V_set;
//         for(int node: temp_set){

//             double dis_1 = euclidean_distance(dataset[min_idx], dataset[node]);
//             auto it = find_in_vector(visited_nodes, node);
//             double dis_2 = it->second;

//             if(alpha * dis_1 <= dis_2){
//                 visited_nodes.erase(it);
//             }
//         }
//     }

//     //update the real graph neighbors with the new N_out neighbors 
//     graph.get_node(query)->neighbors = N_out;

//     return;
// }




// double find_distance_in_map(std::map<double, int> map, int node){

//     for(const auto& pair : map){
//         if(pair.second == node){
//             double key = pair.first; 
//             return key;
//         }
//     }

//     return -1.0;
// }


// template <typename type>
// void Vamana::RobustPruning(RRGraph graph, int query, std::unordered_set<int> V_set, float alpha, int R, std::vector<std::vector<type> > dataset){

//     std::vector<int> N_out = graph.get_node(query)->neighbors;
//     V_set.insert(N_out.begin(), N_out.end());
//     V_set.erase(query);

//     N_out.clear();

//     std::map<double, int> distance_map;

//     for(int node : V_set){
        
//         double dis = euclidean_distance(dataset[query], dataset[node]);
//         distance_map.insert({dis, node});
//     }

//     while(!V_set.empty()){      

//         auto it = distance_map.begin();
//         int min_idx = it->second;

//         //if(std::find(N_out.begin(), N_out.end(), min_idx) == N_out.end()){
//             N_out.push_back(min_idx);
//         //}

//         if(N_out.size() == R) break;

//         std::unordered_set<int> temp_set = V_set;
//         for(int node: temp_set){

//             double dis_1 = euclidean_distance(dataset[min_idx], dataset[node]);
//             double dis_2 = find_distance_in_map(distance_map, node);

//             if(alpha * dis_1 <= dis_2){
//                 V_set.erase(node);
//                 distance_map.erase(dis_2);
//             }
//         }

//     }

//     //update the real graph neighbors with the new N_out neighbors 
//     graph.get_node(query)->neighbors = N_out;

//     return;
// }