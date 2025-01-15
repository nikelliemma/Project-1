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
#include <unordered_map>
#include <future>
#include <mutex> 
#include <omp.h>

#include "Vamana.h"
#include "Graph.h"
#include "dataset.h"
#include "filtered_dataset.h"




//constructor
Vamana::Vamana(int R, int L, int alpha): vamana_index(R){
    this->filename = "";
    this->alpha = alpha;
    this->L = L;
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

void Vamana::set_index(RRGraph graph){
    this->vamana_index = graph;
}

RRGraph Vamana::get_index(){
    return this->vamana_index;
}

double min_distance = 0.0;

//templated function to compute Euclidean distance between two vectors
template <typename type>
double Vamana::euclidean_distance(const std::vector<type>& vec1, const std::vector<type>& vec2){

    if(vec1.size() != vec2.size()){
        throw std::invalid_argument("Vectors must be of the same length.");
    }

    double dist = 0.0;

    //#pragma omp parallel for reduction(+:dist) 
    for(int i = 0; i < vec1.size(); ++i){

        double diff = vec1[i] - vec2[i];
        dist += diff * diff;

        // /if(dist > min_distance) return std::numeric_limits<double>::max();

    }

    return std::sqrt(dist);
}



double euclidean_distance(std::vector<float>& vec1, std::vector<float>& vec2){

    if(vec1.size() != vec2.size()){
        throw std::invalid_argument("Vectors must be of the same length.");
    }

    double dist = 0.0;
    
    //#pragma omp parallel for reduction(+:dist) 
    for(int i = 0; i < vec1.size(); ++i){

        double diff = vec1[i] - vec2[i];
        dist += diff * diff;

        if(dist > min_distance) return MAX_DIS;

    }

    return std::sqrt(dist);
}

double euclidean_distance1(std::vector<float>& vec1, std::vector<float>& vec2){

    if(vec1.size() != vec2.size()){
        throw std::invalid_argument("Vectors must be of the same length.");
    }

    double dist = 0.0;
    
    //#pragma omp parallel for reduction(+:dist) 
    for(int i = 0; i < vec1.size(); ++i){

        double diff = vec1[i] - vec2[i];
        dist += diff * diff;

    }

    return std::sqrt(dist);
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

int Vamana::find_medoid_filtered(std::vector<Data_Point> dataset, int filter) {

    int n = dataset.size();
    int medoid_index = -1;
    double min_total_distance = std::numeric_limits<double>::max();

    // Filtered indices of the dataset that match the categorical value
    std::vector<int> filtered_indices;
    for (int i = 0; i < n; ++i) {
        if (dataset[i].categorical == filter) {
            filtered_indices.push_back(i);
        }
    }

    int filtered_size = filtered_indices.size();

    // Return -1 if no data points match the filter
    if (filtered_size == 0) {
        return -1;
    }

    // Compute the medoid using only filtered data points
    for (int i_idx = 0; i_idx < filtered_size; ++i_idx) {
        int i = filtered_indices[i_idx];
        double total_distance = 0.0;

        for (int j_idx = 0; j_idx < filtered_size; ++j_idx) {
            int j = filtered_indices[j_idx];
            if (i != j) {
                total_distance += euclidean_distance(dataset[i].data_vector, dataset[j].data_vector);
            }
        }

        if (total_distance < min_total_distance) {
            min_total_distance = total_distance;
            medoid_index = i;
        }
    }

    return medoid_index;
}





int Vamana::find_medoid_f(std::vector<Data_Point> dataset){

    int n = dataset.size();
    int d = dataset[0].data_vector.size();
    int medoid_index = -1;
    double min_total_distance = std::numeric_limits<double>::max();
    
    for(int i = 0; i < n; ++i){

        double total_distance = 0.0;
        
        for(int j = 0; j < n; ++j){
            if(i != j){
                total_distance += euclidean_distance(dataset[i].data_vector, dataset[j].data_vector);
            }
        }
        
        if(total_distance < min_total_distance){
            min_total_distance = total_distance;
            medoid_index = i;
        }
    }
    
    return medoid_index;    
}


//check if the min heap has elements that are not in the visited set
bool has_unvisited_elements(std::priority_queue<std::pair<double, int>, 
                        std::vector<std::pair<double, int>>, 
                        std::greater<std::pair<double, int>>>& min_heap, 
                        const std::unordered_set<int>& visited){
    
    //create a copy of the min_heap to avoid modifying the original heap
    auto min_heap_copy = min_heap;

    //iterate over the elements in the min_heap
    while(!min_heap_copy.empty()){

        //get the top element (minimum element in the min-heap)
        int node_id = min_heap_copy.top().second; 
        min_heap_copy.pop(); 

        //check if this element is not in the visited set
        if(visited.find(node_id) == visited.end()){
            return true; //found an unvisited elementv
        }
    }

    return false; //all elements in the min_heap are visited
}


//greedy search function
template <typename type>
LVPair Vamana::GreedySearch(RRGraph graph, int starting_node, int query, int k, int L, std::vector<std::vector<type> > dataset){

    //result vector contains the NNs of the query
    //visited set contains the visited nodes 
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
    std::unordered_set<int> visited;
    std::vector<int> result;

    //set to keep track of the nodes that are in min_heap
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
            //min_heap_nodes.clear();
            min_heap = {};
            for(int i = 0; i < L; i++){
                min_heap.emplace(min_heap_copy.top());
                //min_heap_nodes.emplace(min_heap_copy.top().second);
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

//second greedy search function that instead of the index of the node it takes the whole vector as an argument 
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

            //insert them only if they are not already in the min heap 
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
            //min_heap_nodes.clear();
            min_heap = {};
            for(int i = 0; i < L; i++){
                min_heap.emplace(min_heap_copy.top());
                //min_heap_nodes.emplace(min_heap_copy.top().second);
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
    
    //return the k-NNs and the visited nodes set
    return {result, visited};
}






bool has_unvisited_elements(std::vector<std::pair<int,double>> vec, std::unordered_set<int> set){

    for(auto node: vec){
        //check if this element is not in the visited set
        if(set.find(node.first) == set.end()){
            return true; //found an unvisited element
        }
    }

    return false;

}

LVPair Vamana::FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset){

    //result vector contains the NNs of the query
    //visited set contains the visited nodes 
    std::vector<std::pair<int,double>> L_set;
    std::unordered_set<int> visited;
    std::vector<int> result;

    for(const auto & [filter, node]: S_nodes){
        if(filter == dataset[query_vec].categorical){
            L_set.push_back({node, euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector)});
        }
    }

    //start the greedy search loop
    while(has_unvisited_elements(L_set, visited)){

        //get the minimum distance of the nodes in L set that have node been visited

        std::sort(L_set.begin(), L_set.end(),
            [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second < b.second; // Ascending order
            });

        int index = 0;
        for(auto node : L_set){
            if(visited.find(node.first) == visited.end()) break;
            index++;
        }
        int min_idx = L_set[index].first;
        min_distance = L_set[index].second;

        std::vector<int> new_neighbours;

        //update the visited set of visited nodes
        visited.insert(min_idx);
       
        //insert the neighbors of the min distance node to the query into the result (L) set
        for(int node: graph.get_node(min_idx)->neighbors){
   
            Data_Point point = dataset[node];

            if(point.categorical == dataset[min_idx].categorical){
                if(visited.find(node) == visited.end()){
                    new_neighbours.push_back(node);
                }
            }

        }


        for(int node: new_neighbours){

            auto it = std::find_if(L_set.begin(), L_set.end(), [node](const std::pair<int, double>& p) {
                return p.first == node;
            });

            // If the node is not found in L_set, add it
            if(it == L_set.end()){
                double dist = euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector);
                if (dist < min_distance) {
                    min_distance = dist; // Update the minimum distance
                }
                // L_set.push_back({node, euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector)});
                L_set.push_back({node, dist});
            }
        }

        //retain the L closest points to the query
        if(L_set.size() > L){

            std::sort(L_set.begin(), L_set.end(),
                [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                    return a.second < b.second; // Ascending order
                });

            L_set.resize(L);

        }

    }

    std::sort(L_set.begin(), L_set.end(),[](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return a.second < b.second; // Ascending order
    });
    
    int counter = 0;

    //get the k nearest neighbors
    for(auto node: L_set){
        result.push_back(node.first);
        if(counter == k) break;
    }
    
    //return the k-NNs and the visited nodes set
    return {result, visited};
    
}


LVPair Vamana::FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset, std::vector<Data_Point> Q_dataset){

    //result vector contains the NNs of the query
    //visited set contains the visited nodes 
    std::vector<std::pair<int,double>> L_set;
    std::unordered_set<int> visited;
    std::vector<int> result;

    for(const auto & [filter, node]: S_nodes){
        if(filter == Q_dataset[query_vec].categorical || Q_dataset[query_vec].categorical == -1){
            L_set.push_back({node, euclidean_distance1(dataset[node].data_vector, Q_dataset[query_vec].data_vector)});
        }
    }

    //start the greedy search loop
    while(has_unvisited_elements(L_set, visited)){

        //get the minimum distance of the nodes in L set that have node been visited

        std::sort(L_set.begin(), L_set.end(),
            [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second < b.second; // Ascending order
            });

        int index = 0;
        for(auto node : L_set){
            if(visited.find(node.first) == visited.end()) break;
            index++;
        }
        int min_idx = L_set[index].first;
        min_distance = L_set[index].second;

        std::vector<int> new_neighbours;

        //update the visited set of visited nodes
        visited.insert(min_idx);
       
        //insert the neighbors of the min distance node to the query into the result (L) set
        for(int node: graph.get_node(min_idx)->neighbors){
   
            Data_Point point = dataset[node];

            if(point.categorical == dataset[min_idx].categorical){
                if(visited.find(node) == visited.end()){
                    new_neighbours.push_back(node);
                }
            }

        }

        for(int node: new_neighbours){

            auto it = std::find_if(L_set.begin(), L_set.end(), [node](const std::pair<int, double>& p) {
                return p.first == node;
            });

            // If the node is not found in L_set, add it
            if (it == L_set.end()) {
                L_set.push_back({node, euclidean_distance1(dataset[node].data_vector, Q_dataset[query_vec].data_vector)});
            }
        }

        //retain the L closest points to the query
        if(L_set.size() > L){

            std::sort(L_set.begin(), L_set.end(),
                [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                    return a.second < b.second; // Ascending order
                });

            L_set.resize(L);

        }

    }

    std::sort(L_set.begin(), L_set.end(),[](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return a.second < b.second; // Ascending order
    });
    
    int counter = 0;

    //get the k nearest neighbors
    for(auto node: L_set){
        result.push_back(node.first);
        if(counter == k) break;
    }
    
    //return the k-NNs and the visited nodes set
    return {result, visited};
    
}






//function that generates a random permutation from 0 to n-1
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


std::vector<int> get_ids(std::vector<Data_Point> dataset, int filter){

    std::vector<int> P_set;

    int size = dataset.size();

    for(int i = 0; i<size ;i++){
        int fil = dataset[i].categorical;

        if(fil == filter || fil == -1){
            P_set.push_back(i);
        }

    }
    
    return P_set;
}

std::vector<int> select_random_elements(const std::vector<int> original, int t){

    // if(t > original.size()){
    //     //throw std::invalid_argument("t cannot be greater than the size of the original vector.");
    // }

    std::vector<int> shuffled = original;
    std::random_device rd; // Obtain a random number generator
    std::mt19937 gen(rd()); // Seed the generator
    std::shuffle(shuffled.begin(), shuffled.end(), gen);

    if(t > original.size()) return std::vector<int>(shuffled.begin(), shuffled.begin() + 1);
    return std::vector<int>(shuffled.begin(), shuffled.begin() + t);

}


//parallel find filtered medoid
std::map<int, int> Vamana::Filtered_Find_Medoid(std::vector<Data_Point> dataset, std::unordered_set<int> filters, int threshold) {
    std::map<int, int> M_map; 
    std::map<int, int> T_map; 

    //#pragma omp parallel for schedule(runtime)
    for(int i = 0; i < dataset.size(); ++i){
        T_map[i] = 0;
    }

    //iterate over each filter
    #pragma omp parallel for schedule(runtime)
    for(int filter = 0; filter < filters.size(); filter++){
    //for(const int &filter : filters){
    
        std::vector<int> F_ids = get_ids(dataset, filter);

        std::vector<int> random_points = select_random_elements(F_ids, threshold);

        //find the point with the minimum usage in T_map
        int min_count = std::numeric_limits<int>::max();
        int min_point = -1;

        for(int point : random_points){  
            #pragma omp critical
            {
                if(T_map[point] < min_count){
                    min_count = T_map[point];
                    min_point = point;
                }
            }
        }

        
        //update M_map and T_map
        if(min_point != -1){ 
            #pragma omp critical
            {
                M_map[filter] = min_point;
                T_map[min_point] += 1;
            }
        }
    }

    return M_map;
}



//robust pruning algorithm
template <typename type>
void Vamana::RobustPruning(RRGraph G, int q, std::unordered_set<int> V, float a, int R, std::vector<std::vector<type>> dataset){
    std::vector<int> out = G.get_node(q)->neighbors;
    V.insert(out.begin(), out.end());
    V.erase(q);
    out.clear();

    //create minheap of nodes and their distance from q
    std::priority_queue<std::pair<double, int>> minHeap{};

    for(int node : V){
        minHeap.emplace(-1 * euclidean_distance(dataset[q], dataset[node]), node);//* (-1) to make it min from max
    }

    while(!V.empty()){
        int minNode = minHeap.top().second;//get closest node to q

        if(std::find(out.begin(), out.end(), minNode) == out.end()){ //if not there add it to out
            out.push_back(minNode);
        }

        //if we have enough nodes, stop
        if(out.size() == R){
            break;
        }

        std::unordered_set<int> temp = V;
        
        for(int node: temp){
            double dis1 = euclidean_distance(dataset[minNode], dataset[node]);
            double dis2 = euclidean_distance(dataset[q], dataset[node]);


            if((a * dis1) <= dis2){
                V.erase(node);

                //remake heap for updated V
                minHeap = {};
                for(int i : V){
                    minHeap.emplace(-1 * euclidean_distance(dataset[q], dataset[i]), i);
                }
            }
        }

    }

    //update new neighboors
    G.get_node(q)->neighbors = std::move(out);

}


// filtered robust pruning algorithm
void Vamana::FilteredRobustPruning(RRGraph G, int q, std::unordered_set<int> V, float a, int R, FilteredDataset f_dataset){

    std::vector<Data_Point> dataset = f_dataset.get_dataset();

    std::vector<int> out = G.get_node(q)->neighbors;
    V.insert(out.begin(), out.end());
    V.erase(q);
    out.clear();

    //create minheap of nodes and their distance from q
    std::priority_queue<std::pair<double, int>> minHeap{};

    for(int node : V){
        minHeap.emplace(-1 * euclidean_distance1(dataset[q].data_vector, dataset[node].data_vector), node);//* (-1) to make it min from max
    }

    while(!V.empty()){
       int minNode = minHeap.top().second;//get closest node to q

       if(std::find(out.begin(), out.end(), minNode) == out.end()){ //if not there add it to out
            out.push_back(minNode);
        }

        //if we have enough nodes, stop
        if(out.size() == R){
            break;
        }

        std::unordered_set<int> temp = V;
        

        for(int node: temp){
            
            int filterNode = dataset[node].categorical;
            int filterQ = dataset[q].categorical;
            int filterMin = dataset[minNode].categorical;
        
            if(filterNode == filterQ && filterNode != filterMin){
                continue;
            }
    
            double dis1 = euclidean_distance1(dataset[minNode].data_vector, dataset[node].data_vector);
            double dis2 = euclidean_distance1(dataset[q].data_vector, dataset[node].data_vector);

            if((a * dis1) <= dis2){
                V.erase(node);

                //remake heap for updated V
                minHeap = {};
                for(int i : V){
                   minHeap.emplace(-1 * euclidean_distance1(dataset[q].data_vector, dataset[i].data_vector), i);
                }
            }
        }

    }

    //update new neighboors
    G.get_node(q)->neighbors = std::move(out);
}


// PARALLEL 1

// RRGraph Vamana::Filtered_Vamana_Index_Parallel(FilteredDataset dataset_obj, int L, int R, float a) {
//     std::unordered_set<int> filters_set = dataset_obj.get_filter_set();
//     std::vector<Data_Point> dataset = dataset_obj.get_dataset();

//     //int max_threads = omp_get_num_procs(); // Number of logical cores
//     //omp_set_num_threads(max_threads);     // Set thread count accordingly
//     omp_set_num_threads(4); // Set number of threads
   
    
//     // Create an empty graph
//     RRGraph graph(R);
//     graph.create_Rregular_empty_graph(dataset);
    
//     int N = dataset.size();
    
//     // Get the random permutation as a starting
//     std::vector<int> perm = std::move(get_random_permutation(N)); 

//     // Filter map
//     std::map<int, int> filter_map = Filtered_Find_Medoid(dataset, filters_set, 1);
    
//     // Initialize the visited set
//     std::vector<int> visited;
    
//     // Use OpenMP to parallelize the loop
//     #pragma omp parallel 
//     {
//         // Local private variables for each thread
//         std::vector<int> local_visited;
        
//         #pragma omp for //schedule(dynamic)
//         for (int i = 0; i < N; i++) {
//             cout << i << endl;
//             // Perform greedy search
//             LVPair greedy_result = FilteredGreedySearch(graph, filter_map, perm[i], 0, L, filters_set, dataset);
            
//             // Insert the result into the local_visited list
//             local_visited.insert(local_visited.end(), greedy_result.second.begin(), greedy_result.second.end());

//             // Run filtered robust prune
//             FilteredRobustPruning(graph, perm[i], greedy_result.second, a, R, dataset_obj);

//             // Get the neighbors of perm[i]
//             std::vector<int> perm_i_neighbours = graph.get_node(perm[i])->neighbors;

//             for (int j : perm_i_neighbours) {
//                 // Ensure mutual neighbors (i.e., undirected graph behavior)
//                 if (std::find(graph.get_node(j)->neighbors.begin(), graph.get_node(j)->neighbors.end(), perm[i]) == graph.get_node(j)->neighbors.end()) {
//                     #pragma omp critical
//                     graph.get_node(j)->neighbors.push_back(perm[i]);
//                 }

//                 int n = graph.get_node(j)->neighbors.size();
//                 std::unordered_set<int> N_out_j(graph.get_node(j)->neighbors.begin(), graph.get_node(j)->neighbors.end());

//                 if (n > R) {
//                     // Run filtered robust prune on the neighbors if needed
//                     FilteredRobustPruning(graph, j, N_out_j, a, R, dataset_obj);
//                 }
//             }
//         }
        
//         // Now, merge local_visited back to the global visited list safely
//         #pragma omp critical
//         {
//             visited.insert(visited.end(), local_visited.begin(), local_visited.end());
//         }
//     }

//     return graph;
// }


// PARALLEL 2

RRGraph Vamana::Filtered_Vamana_Index_Parallel(FilteredDataset dataset_obj, int L, int R, float a){

    std::unordered_set<int> filters_set = dataset_obj.get_filter_set();
    std::vector<Data_Point> dataset = dataset_obj.get_dataset();

    int max_threads = omp_get_num_procs(); // Number of logical cores
    omp_set_num_threads(max_threads);     // Set thread count accordingly
    //omp_set_num_threads(4); // Set number of threads
    
    // Create an empty graph
    RRGraph graph(R);
    graph.create_Rregular_empty_graph(dataset);
    
    int N = dataset.size();
    
    // Get the random permutation as a starting
    std::vector<int> perm = std::move(get_random_permutation(N)); 

    // Filter map
    std::map<int, int> filter_map = Filtered_Find_Medoid(dataset, filters_set, 1);
    
    // Initialize the visited set
    std::vector<int> visited;
    
    // Use OpenMP to parallelize the loop
    #pragma omp parallel
    {
        // Local private variables for each thread
        std::vector<int> local_visited;
        
        #pragma omp for nowait
        for (int i = 0; i < N; i++) {

            //cout << i << endl;
            // Perform greedy search
            LVPair greedy_result = FilteredGreedySearch(graph, filter_map, perm[i], 0, L, filters_set, dataset);
            
            // Insert the result into the local_visited list
            local_visited.insert(local_visited.end(), greedy_result.second.begin(), greedy_result.second.end());

            // Run filtered robust prune
            FilteredRobustPruning(graph, perm[i], greedy_result.second, a, R, dataset_obj);

            // Get the neighbors of perm[i]
            std::vector<int> perm_i_neighbours = graph.get_node(perm[i])->neighbors;

            for (int j : perm_i_neighbours) {
                
                if (std::find(graph.get_node(j)->neighbors.begin(), graph.get_node(j)->neighbors.end(), perm[i]) == graph.get_node(j)->neighbors.end()) {
                    #pragma omp critical
                    graph.get_node(j)->neighbors.push_back(perm[i]);
                }

                int n = graph.get_node(j)->neighbors.size();
                std::unordered_set<int> N_out_j(graph.get_node(j)->neighbors.begin(), graph.get_node(j)->neighbors.end());

                if (n > R) {
                    // Run filtered robust prune on the neighbors if needed
                    FilteredRobustPruning(graph, j, N_out_j, a, R, dataset_obj);
                }
            }
        }
        
        #pragma omp critical
        {
            visited.insert(visited.end(), local_visited.begin(), local_visited.end());
        }
    }

    return graph;
}






// PARALLEL 3 

// RRGraph Vamana::Filtered_Vamana_Index_Parallel(FilteredDataset dataset_obj, int L, int R, float a){

//     std::unordered_set<int> filters_set = dataset_obj.get_filter_set();
//     std::vector<Data_Point> dataset = dataset_obj.get_dataset();

//     omp_set_num_threads(4);
    
//     //create an empty graph
//     RRGraph graph(R);
//     graph.create_Rregular_empty_graph(dataset);
//     //graph.create_Rregular_graph(dataset); it makes it too slow!! 

//     //get the size of the dataset()
//     int N = dataset.size();

//     //get the random permutation as a starting 
//     std::vector<int> perm = std::move(get_random_permutation(N)); 

//     std::map<int, int> filter_map = Filtered_Find_Medoid(dataset, filters_set, 1);

//     //L = filter_map.size() * 10;

//     std::vector<int> visited;

//     #pragma omp parallel
//     {

//     #pragma omp for schedule(dynamic)       
//     for(int i = 0; i<N; i++){

//         cout << i << endl;

//         LVPair greedy_result = FilteredGreedySearch(graph, filter_map, perm[i], 0, L, filters_set, dataset);

//         #pragma omp critical
//         visited.insert(visited.end(), greedy_result.second.begin(), greedy_result.second.end());

//         //run filtered robust prune
//         #pragma omp critical
//         FilteredRobustPruning(graph, perm[i], greedy_result.second, a, R, dataset_obj);

//         //get the neighbors of perm[i]
//         std::vector<int> perm_i_neighbours = graph.get_node(perm[i])->neighbors;

//         for(int j: perm_i_neighbours){
            
//             if(std::find((graph.get_node(j)->neighbors).begin(), (graph.get_node(j)->neighbors).end(), perm[i]) == (graph.get_node(j)->neighbors).end()) {
//                 #pragma omp critical 
//                 (graph.get_node(j)->neighbors).push_back(perm[i]);
//             }

//             int n = graph.get_node(j)->neighbors.size();

//             std::unordered_set<int> N_out_j(graph.get_node(j)->neighbors.begin(), graph.get_node(j)->neighbors.end());

//             if(n > R){
//                 //run filtered robust prune
//                 #pragma omp critical
//                 FilteredRobustPruning(graph, j, N_out_j, a, R, dataset_obj);
//             }

//         }

//     }
//     }

//     return graph;

// }


RRGraph Vamana::Filtered_Vamana_Index(FilteredDataset dataset_obj, int L, int R, float a){

    std::unordered_set<int> filters_set = dataset_obj.get_filter_set();
    std::vector<Data_Point> dataset = dataset_obj.get_dataset();

    omp_set_num_threads(16);
    
    //create an empty graph
    RRGraph graph(R);
    graph.create_Rregular_empty_graph(dataset);
    //graph.create_Rregular_graph(dataset); it makes it too slow!! 

    //get the size of the dataset()
    int N = dataset.size();

    //get the random permutation as a starting 
    std::vector<int> perm = std::move(get_random_permutation(N)); 

    std::map<int, int> filter_map = Filtered_Find_Medoid(dataset, filters_set, 1);

    //L = filter_map.size() * 10;

    std::vector<int> visited;

    for(int i = 0; i<N; i++){

        cout << i << endl;

        LVPair greedy_result = FilteredGreedySearch(graph, filter_map, perm[i], 0, L, filters_set, dataset);

        visited.insert(visited.end(), greedy_result.second.begin(), greedy_result.second.end());

        //run filtered robust prune
        FilteredRobustPruning(graph, perm[i], greedy_result.second, a, R, dataset_obj);

        //get the neighbors of perm[i]
        std::vector<int> perm_i_neighbours = graph.get_node(perm[i])->neighbors;

        for(int j: perm_i_neighbours){
            
            if(std::find((graph.get_node(j)->neighbors).begin(), (graph.get_node(j)->neighbors).end(), perm[i]) == (graph.get_node(j)->neighbors).end()) {
                (graph.get_node(j)->neighbors).push_back(perm[i]);
            }

            int n = graph.get_node(j)->neighbors.size();

            std::unordered_set<int> N_out_j(graph.get_node(j)->neighbors.begin(), graph.get_node(j)->neighbors.end());

            if(n > R){
                //run filtered robust prune
                FilteredRobustPruning(graph, j, N_out_j, a, R, dataset_obj);
            }

        }

    }

    return graph;

}




template <typename type>
RRGraph Vamana::Vamana_Index(std::vector<std::vector<type> > dataset, int L, int R, float a){

    //create a R-regular graph based on the geometric properties of the dataset
    RRGraph graph(R);
    graph.create_Rregular_graph(dataset);

    //get the size of the dataset()
    int N = dataset.size();

    //get the medoid of the dataset 
    int medoid_index = find_medoid(dataset);

    //get the random permutation as a starting 
    std::vector<int> perm = std::move(get_random_permutation(N)); 

    //start the Vamana loop 
    for(int i = 0; i < N; i++){

        //cout << i << endl;

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

    for(int i = 0;i<graph.get_nodes_num();i++){
        if(graph.get_node(i)->neighbors.size() > R) graph.get_node(i)->neighbors.resize(R);
    }

    return graph;

}




std::vector<RRGraph> Vamana::StitchedVamana(FilteredDataset dataset_obj, int Lsmall, int Rsmall, int Rstitched, float a){

    std::vector<RRGraph> filter_graphs;
    int n = dataset_obj.get_filter_set().size();
    std::vector<Data_Point> dataset = dataset_obj.get_dataset();

    for(int i = 0; i < n; i++) {
        RRGraph subgraph(Rsmall);
        filter_graphs.push_back(subgraph);
    }

    for(int filter = 0; filter < n; filter++){

        cout << filter << endl;

        std::vector<std::vector<float>> PointsByFilter;

        for(int i = 0; i < dataset.size(); i++){
            if(dataset_obj.get_data_point(i).categorical == filter) PointsByFilter.push_back(dataset_obj.get_data_point(i).data_vector);
        }

        std::cout << PointsByFilter.size() << std::endl;

        if(PointsByFilter.size() <= Rstitched) continue;

        filter_graphs[filter] = Vamana_Index(PointsByFilter, Lsmall, Rstitched, a);

    }

    return filter_graphs;
}




std::vector<RRGraph> Vamana::StitchedVamanaParallel(FilteredDataset dataset_obj, int Lsmall, int Rsmall, int Rstitched, float a){

    std::vector<RRGraph> filter_graphs;
    int n = dataset_obj.get_filter_set().size();
    std::vector<Data_Point> dataset = dataset_obj.get_dataset();

    for(int i = 0; i < n; i++) {
        RRGraph subgraph(Rsmall);
        filter_graphs.push_back(subgraph);
    }

    //use omp to parallelize the loop
    #pragma omp parallel for
    for(int filter = 0; filter < n; filter++){

        // cout << filter << endl;

        std::vector<std::vector<float>> PointsByFilter;

        for(int i = 0; i < dataset.size(); i++){
            if(dataset_obj.get_data_point(i).categorical == filter) PointsByFilter.push_back(dataset_obj.get_data_point(i).data_vector);
        }

        // std::cout << PointsByFilter.size() << std::endl;

        if(PointsByFilter.size() <= Rstitched) continue;

        filter_graphs[filter] = Vamana_Index(PointsByFilter, Lsmall, Rstitched, a);

    }

    return filter_graphs;
}



//function to get the recall 
double Vamana::Get_Recall(std::vector<int> vec1, std::vector<int> vec2){

    if(vec1.size() == 0 || vec2.size() == 0) return 0.0;

    std::unordered_map<int, int> elementCount;
    for(int num : vec1){
        elementCount[num]++;
    }

    //count common elements 
    int count = 0;
    for(int num : vec2){
        if(elementCount[num] > 0){
            ++count;
            --elementCount[num]; 
        }
    }

    //return recall 
    return (count * 100) / vec1.size();

}


int count_filter_points(FilteredDataset dataset, int filter){

    int points = 0;

    for(int i = 0; i < dataset.get_dataset().size(); i++){
        if(dataset.get_data_point(i).categorical == filter) points++;
    }

    return points;
}

double Vamana::Get_Recall_Filtered(std::vector<int> knns, std::vector<int> vec2, FilteredDataset dataset, int k, int filter){

    if(knns.size() == 0 || vec2.size() == 0) return 0.0;

    std::unordered_map<int, int> elementCount;
    for(int num : knns){
        elementCount[num]++;
    }

    //count common elements 
    int count = 0;
    for(int num : vec2){
        if(elementCount[num] > 0){
            ++count;
            --elementCount[num]; 
        }
    }

    int filter_points = count_filter_points(dataset, filter);

    if(k > filter_points && filter_points != 0){
        return count * 100 / filter_points;
    }

    return count * 100 / k;
}



std::string extract_format(std::string &file){

    std::string format;
    int len = file.length();

    int pos = file.find_last_of("/");
    std::string file_string = file.substr(pos + 1);

    int index = 0;
    while(file[index] != '.'){ index++; }

    if(index == len - 1) return "";

    format = file.substr(index + 1);

    return format;

}


//function to identify the given file and create the Vamana Index 
void Vamana::create_vamana_index(std::string filepath, int L, int R, int alpha){

    std::string file_format = extract_format(filepath);

    if(file_format == "fvecs"){
        Dataset<float> d;
        d.set_filename(filepath);
        d.read_dataset();

        if(d.get_dataset().empty()){
            return;
        }
        RRGraph Vam = Vamana_Index(d.get_dataset(), L, R, alpha);
        this->vamana_index = Vam;
        this->vamana_index.print_graph();
    }
    else if(file_format == "ivecs"){
        Dataset<int> d;
        d.set_filename(filepath);
        d.read_dataset();

        if(d.get_dataset().empty()){
            return;
        }
        RRGraph Vam = Vamana_Index(d.get_dataset(), L, R, alpha);   
        this->vamana_index = Vam;
        this->vamana_index.print_graph();
    }
    else if(file_format == "bvecs"){
        Dataset<unsigned char> d;
        d.set_filename(filepath);
        d.read_dataset();

        if(d.get_dataset().empty()){
            return;
        }
        RRGraph Vam = Vamana_Index(d.get_dataset(), L, R, alpha);
        this->vamana_index = Vam;
        this->vamana_index.print_graph();
    }

    return;
}



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













//THIS IS SOME CODE OF SOME OPTIMIZATIONS THAT DIDN'T WORK 


// GraphCollection Vamana::StitchedVamana(FilteredDataset dataset_obj, int Lsmall, int Rsmall, int Rstitched, float a){
    
//     GraphCollection collectionOfGraphs;
//     std::vector<std::vector<float>> PointsByFilter;

//     std::unordered_set<int> filters = dataset_obj.get_filter_set();

//     int size = dataset_obj.get_dataset().size();
//     int fsize = filters.size();

//     RRGraph subGraph(Rsmall);

//     //seperate points by filter and create subgraphs for each
//     for(int f = 0; f < fsize; f++){

//         cout << f << endl;

//         PointsByFilter.clear();

//         //create vector of points for each filter        
//         for(int i = 0; i < size; i++){

//             if(dataset_obj.get_data_point(i).categorical == f){
//                 PointsByFilter.push_back(dataset_obj.get_data_point(i).data_vector);
//             }
//         }
//         //cout << "nodes num : " << PointsByFilter.size() << endl;

//         //create subgraphs with vamana for each filter
//         //since there is no stitching, we run vamana with Rstitched instead of Rsmall
//         if(PointsByFilter.size() <= Rsmall){
//             //skipping 1 or less point graph cause no point?
//             continue;
//         }
        
//         subGraph = Vamana_Index(PointsByFilter, Lsmall, Rsmall, a);
//         collectionOfGraphs.push_back(subGraph);
//     }

    
//     // for (int g = 0; g < collectionOfGraphs.size(); g++){
//     //     int sizeG = collectionOfGraphs[g].get_graph().size();
//     //     for(int q = 0; q < sizeG; q++){
//     //         FilteredRobustPruning(collectionOfGraphs[g], q, , a, Rstitched, dataset_obj);
//     //     }
//     // }
    

//     return collectionOfGraphs;
// }




// LVPair Vamana::FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset){

//     //result vector contains the NNs of the query
//     //visited set contains the visited nodes 
//     std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
//     std::unordered_set<int> visited;
//     std::vector<int> result;
//     std::unordered_set<int> min_heap_nodes;

//     for(const auto & [filter, node]: S_nodes){
//         if(filter == dataset[query_vec].categorical){
//             min_heap.emplace(euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector), node);  
//             min_heap_nodes.emplace(node); 
//         }
//     }

//     //start the greedy search loop
//     while(has_unvisited_elements(min_heap, visited)){

//         //get the minimum distance of the nodes in L set that have node been visited
//         auto min_heap_copy = min_heap;
//         int temp_node;

//         while(!min_heap_copy.empty()){

//             temp_node = min_heap_copy.top().second;
//             if(visited.find(temp_node) != visited.end()){
//                 min_heap_copy.pop();
//             }
//             else{ break; }
//         }

//         int min_idx = temp_node;

//         //update the visited set of visited nodes
//         visited.insert(min_idx);

//         std::vector<int> new_neighbours;
       
//         //insert the neighbors of the min distance node to the query into the result (L) set
//         for(int node: graph.get_node(min_idx)->neighbors){

//             Data_Point point = dataset[node];

//             if(point.categorical == dataset[min_idx].categorical){
//                 if(visited.find(node) == visited.end()){
//                     new_neighbours.push_back(node);
//                 }
//             }

//         }


//         for(int node: new_neighbours){
//             if(min_heap_nodes.find(node) == min_heap_nodes.end()){ 
//                 double distance = euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector);
//                 min_heap.emplace(distance, node);
//                 min_heap_nodes.emplace(node);
//             }
//         }

//         //update the visited set of visited nodes
//         visited.insert(min_idx);

//         //retain the L closest points to the query
//         if(min_heap.size() > L){
//             auto min_heap_copy = min_heap;
//             min_heap = {};
//             min_heap_nodes.clear();
//             for(int i = 0; i < L; i++){
//                 min_heap.emplace(min_heap_copy.top());
//                 min_heap_nodes.emplace(min_heap_copy.top().second);
//                 min_heap_copy.pop();
//             }    
//         }

//     }

//     //get the k nearest neighbors
//     int i = 0;
//     while(i < k){
//         result.push_back(min_heap.top().second);
//         min_heap.pop();
//         i++;
//     }
    
//     //return the k-NNs and the visited nodes set
//     return {result, visited};
    
// }



// LVPair Vamana::FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset, std::vector<Data_Point> Q_dataset){

//     //result vector contains the NNs of the query
//     //visited set contains the visited nodes 
//     std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
//     std::unordered_set<int> visited;
//     std::vector<int> result;
//     std::unordered_set<int> min_heap_nodes;

//     for(const auto & [filter, node]: S_nodes){
//         //if(filter == Q_dataset[query_vec].categorical || Q_dataset[query_vec].categorical == -1){
//         if(filter == Q_dataset[query_vec].categorical){
//             min_heap.emplace(euclidean_distance(dataset[node].data_vector, Q_dataset[query_vec].data_vector), node);  
//             min_heap_nodes.emplace(node);  
//         }
//     }

//     //start the greedy search loop
//     while(has_unvisited_elements(min_heap, visited)){

//         //get the minimum distance of the nodes in L set that have node been visited
//         auto min_heap_copy = min_heap;
//         int temp_node;

//         while(!min_heap_copy.empty()){

//             temp_node = min_heap_copy.top().second;
//             if(visited.find(temp_node) != visited.end()){
//                 min_heap_copy.pop();
//             }
//             else{ break; }
//         }

//         int min_idx = temp_node;

//         //update the visited set of visited nodes
//         visited.insert(min_idx);

//         std::vector<int> new_neighbours;
       
//         //insert the neighbors of the min distance node to the query into the result (L) set
//         for(int node: graph.get_node(min_idx)->neighbors){

//             Data_Point point = dataset[node];

//             //if(point.categorical == dataset[min_idx].categorical || point.categorical == -1){
//             if(point.categorical == dataset[min_idx].categorical){
//                 if(visited.find(node) == visited.end()){
//                     new_neighbours.push_back(node);
//                 }
//             }

//         }


//         for(int node: new_neighbours){
//             if(min_heap_nodes.find(node) == min_heap_nodes.end()){
//                 double distance = euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector);
//                 min_heap.emplace(distance, node);
//                 min_heap_nodes.emplace(node);
//             }
//         }

//         //update the visited set of visited nodes
//         visited.insert(min_idx);

//         //retain the L closest points to the query
//         if(min_heap.size() > L){
//             auto min_heap_copy = min_heap;
//             min_heap = {};
//             min_heap_nodes.clear();
//             for(int i = 0; i < L; i++){
//                 min_heap.emplace(min_heap_copy.top());
//                 min_heap_nodes.emplace(min_heap_copy.top().second);
//                 min_heap_copy.pop();
//             }    
//         }

//     }

//     //get the k nearest neighbors
//     int i = 0;
//     while(i < k){
//         result.push_back(min_heap.top().second);
//         min_heap.pop();
//         i++;
//     }
    
//     //return the k-NNs and the visited nodes set
//     return {result, visited};
    
// }





// std::map<int, int> Vamana::Filtered_Find_Medoid(std::vector<Data_Point> dataset, std::unordered_set<int> filters, int threshold) {
//     std::map<int, int> M_map; 
//     std::map<int, int> T_map; 

//     for(int i = 0; i < dataset.size(); ++i) {
//         T_map[i] = 0;
//     }

//     // Pre-compute a global medoid for fallback
//     int global_medoid = 0; // Use a heuristic to find the global medoid (e.g., smallest index or most central point)

//     // Iterate over each filter
//     for(const int &filter : filters) {
//         std::vector<int> F_ids = get_ids(dataset, filter);

//         if(F_ids.empty()) {
//             // Assign global medoid if no points exist for the filter
//             M_map[filter] = global_medoid;
//             continue;
//         }

//         std::vector<int> random_points = select_random_elements(F_ids, threshold);

//         // Find the point with the minimum normalized usage count in T_map
//         double min_score = std::numeric_limits<double>::max();
//         int min_point = -1;

//         for(int point : random_points) {
//             double score = static_cast<double>(T_map[point]) / F_ids.size();
//             if(score < min_score) {
//                 min_score = score;
//                 min_point = point;
//             }
//         }

//         // Update M_map and T_map
//         if(min_point != -1) { 
//             M_map[filter] = min_point;
//             T_map[min_point] += 1;
//         }
//     }

//     return M_map;
// }






// LVPair Vamana::FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset){

//     //result vector contains the NNs of the query
//     //visited set contains the visited nodes 
//     std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
//     std::unordered_set<int> visited;
//     std::vector<int> result;
//     std::unordered_set<int> min_heap_nodes;

//     for(const auto & [filter, node]: S_nodes){
//         if(filter == dataset[query_vec].categorical){
//             min_heap.emplace(euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector), node);  
//             min_heap_nodes.emplace(node); 
//         }
//     }

//     //start the greedy search loop
//     while(has_unvisited_elements(min_heap, visited)){

//         //get the minimum distance of the nodes in L set that have node been visited
//         auto min_heap_copy = min_heap;
//         int temp_node;

//         while(!min_heap_copy.empty()){

//             temp_node = min_heap_copy.top().second;
//             if(visited.find(temp_node) != visited.end()){
//                 min_heap_copy.pop();
//             }
//             else{ break; }
//         }

//         int min_idx = temp_node;

//         //update the visited set of visited nodes
//         visited.insert(min_idx);

//         std::vector<int> new_neighbours;
       
//         //insert the neighbors of the min distance node to the query into the result (L) set
//         for(int node: graph.get_node(min_idx)->neighbors){

//             Data_Point point = dataset[node];

//             if(point.categorical == dataset[min_idx].categorical){
//                 if(visited.find(node) == visited.end()){
//                     new_neighbours.push_back(node);
//                 }
//             }

//         }

//         // //update the visited set of visited nodes
//         // visited.insert(min_idx);

//         for(int node: new_neighbours){
//             if(min_heap_nodes.find(node) == min_heap_nodes.end()){ 
//                 double distance = euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector);
//                 min_heap.emplace(distance, node);
//                 min_heap_nodes.emplace(node);
//             }
//         }

//         //retain the L closest points to the query
//         if(min_heap.size() > L){
//             auto min_heap_copy = min_heap;
//             min_heap = {};
//             min_heap_nodes.clear();
//             for(int i = 0; i < L; i++){
//                 min_heap.emplace(min_heap_copy.top());
//                 min_heap_nodes.emplace(min_heap_copy.top().second);
//                 min_heap_copy.pop();
//             }    
//         }

//     }

//     //get the k nearest neighbors
//     int i = 0;
//     while(i < k){
//         result.push_back(min_heap.top().second);
//         min_heap.pop();
//         i++;
//     }
    
//     //return the k-NNs and the visited nodes set
//     return {result, visited};
    
// }



// LVPair Vamana::FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset, std::vector<Data_Point> Q_dataset){

//     //result vector contains the NNs of the query
//     //visited set contains the visited nodes 
//     std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
//     std::unordered_set<int> visited;
//     std::vector<int> result;
//     std::unordered_set<int> min_heap_nodes;

//     for(const auto & [filter, node]: S_nodes){
//         //if(filter == Q_dataset[query_vec].categorical || Q_dataset[query_vec].categorical == -1){
//         if(filter == Q_dataset[query_vec].categorical){
//             min_heap.emplace(euclidean_distance(dataset[node].data_vector, Q_dataset[query_vec].data_vector), node);  
//             min_heap_nodes.emplace(node);  
//         }
//     }

//     //start the greedy search loop
//     while(has_unvisited_elements(min_heap, visited)){

//         //get the minimum distance of the nodes in L set that have node been visited
//         auto min_heap_copy = min_heap;
//         int temp_node;

//         while(!min_heap_copy.empty()){

//             temp_node = min_heap_copy.top().second;
//             if(visited.find(temp_node) != visited.end()){
//                 min_heap_copy.pop();
//             }
//             else{ break; }
//         }

//         int min_idx = temp_node;

//         //update the visited set of visited nodes
//         visited.insert(min_idx);

//         std::vector<int> new_neighbours;
       
//         //insert the neighbors of the min distance node to the query into the result (L) set
//         for(int node: graph.get_node(min_idx)->neighbors){

//             Data_Point point = dataset[node];

//             //if(point.categorical == dataset[min_idx].categorical || point.categorical == -1){
//             if(point.categorical == dataset[min_idx].categorical){
//                 if(visited.find(node) == visited.end()){
//                     new_neighbours.push_back(node);
//                 }
//             }

//         }

//         // //update the visited set of visited nodes
//         // visited.insert(min_idx);

//         for(int node: new_neighbours){
//             if(min_heap_nodes.find(node) == min_heap_nodes.end()){
//                 double distance = euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector);
//                 min_heap.emplace(distance, node);
//                 min_heap_nodes.emplace(node);
//             }
//         }

//         //retain the L closest points to the query
//         if(min_heap.size() > L){
//             auto min_heap_copy = min_heap;
//             min_heap = {};
//             min_heap_nodes.clear();
//             for(int i = 0; i < L; i++){
//                 min_heap.emplace(min_heap_copy.top());
//                 min_heap_nodes.emplace(min_heap_copy.top().second);
//                 min_heap_copy.pop();
//             }    
//         }

//     }

//     // if(min_heap.size() == 0) return {result, visited};

//     //get the k nearest neighbors
//     int i = 0;
//     while(i < k){
//         result.push_back(min_heap.top().second);
//         min_heap.pop();
//         i++;
//     }
    
//     //return the k-NNs and the visited nodes set
//     return {result, visited};
    
// }