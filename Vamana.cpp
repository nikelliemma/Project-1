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



//templated function to compute Euclidean distance between two vectors
template <typename type>
double Vamana::euclidean_distance(const std::vector<type>& vec1, const std::vector<type>& vec2){

    if(vec1.size() != vec2.size()){
        throw std::invalid_argument("Vectors must be of the same length.");
    }

    double dist = 0.0;
    for(int i = 0; i < vec1.size(); ++i){

        double diff = vec1[i] - vec2[i];
        dist += diff * diff;

    }

    return std::sqrt(dist);
}

double euclidean_distance(std::vector<float>& vec1, std::vector<float>& vec2){

    if(vec1.size() != vec2.size()){
        throw std::invalid_argument("Vectors must be of the same length.");
    }

    double dist = 0.0;
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
            return true; //found an unvisited element
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
    
    //return the k-NNs and the visited nodes set
    return {result, visited};
}


LVPair Vamana::FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset){
//LVPair Vamana::FilteredGreedySearch(RRGraph graph, int starting_node, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset){
    //result vector contains the NNs of the query
    //visited set contains the visited nodes 
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
    std::unordered_set<int> visited;
    std::vector<int> result;
    std::unordered_set<int> min_heap_nodes;

    for(const auto & [filter, node]: S_nodes){
        //if(filter == dataset[query_vec].categorical || filter == -1){
        if(filter == dataset[query_vec].categorical){
            min_heap.emplace(euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector), node);   
        }
    }
    // if(dataset[starting_node].categorical == dataset[query_vec].categorical || dataset[starting_node].categorical == -1 || dataset[query_vec].categorical == -1){
    //     min_heap.emplace(euclidean_distance(dataset[starting_node].data_vector, dataset[query_vec].data_vector), starting_node);
    // }

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

        //update the visited set of visited nodes
        visited.insert(min_idx);

        std::vector<int> new_neighbours;
       
        //insert the neighbors of the min distance node to the query into the result (L) set
        for(int node: graph.get_node(min_idx)->neighbors){

            Data_Point point = dataset[node];

            //if(point.categorical == dataset[min_idx].categorical || point.categorical == -1){
            if(point.categorical == dataset[min_idx].categorical){
                if(visited.find(node) == visited.end()){
                    new_neighbours.push_back(node);
                }
            }

        }

        for(int node: new_neighbours){
            if(min_heap_nodes.find(node) == min_heap_nodes.end()){ 
                double distance = euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector);
                min_heap.emplace(distance, node);
                min_heap_nodes.emplace(node);
            }
        }

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

    //if(min_heap.size() == 0) return {result, visited};

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



LVPair Vamana::FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset, std::vector<Data_Point> Q_dataset){
    //result vector contains the NNs of the query
    //visited set contains the visited nodes 
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> min_heap;
    std::unordered_set<int> visited;
    std::vector<int> result;
    std::unordered_set<int> min_heap_nodes;

    for(const auto & [filter, node]: S_nodes){
        if(filter == Q_dataset[query_vec].categorical || Q_dataset[query_vec].categorical == -1){
        //if(filter == Q_dataset[query_vec].categorical){
            min_heap.emplace(euclidean_distance(dataset[node].data_vector, Q_dataset[query_vec].data_vector), node);   
        }
    }

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

        //update the visited set of visited nodes
        visited.insert(min_idx);

        std::vector<int> new_neighbours;
       
        //insert the neighbors of the min distance node to the query into the result (L) set
        for(int node: graph.get_node(min_idx)->neighbors){

            Data_Point point = dataset[node];

            if(point.categorical == dataset[min_idx].categorical || point.categorical == -1){
            //if(point.categorical == dataset[min_idx].categorical){
                if(visited.find(node) == visited.end()){
                    new_neighbours.push_back(node);
                }
            }

        }

        for(int node: new_neighbours){
            if(min_heap_nodes.find(node) == min_heap_nodes.end()){
                double distance = euclidean_distance(dataset[node].data_vector, dataset[query_vec].data_vector);
                min_heap.emplace(distance, node);
                min_heap_nodes.emplace(node);
            }
        }

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

    //if(min_heap.size() == 0) return {result, visited};

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

        if( fil == -1 || fil == filter ){
            P_set.push_back(i);
        }

    }
    
    return P_set;
}

std::vector<int> select_random_elements(const std::vector<int> original, int t){

    if(t > original.size()){
        throw std::invalid_argument("t cannot be greater than the size of the original vector.");
    }

    std::vector<int> shuffled = original;
    std::random_device rd; // Obtain a random number generator
    std::mt19937 gen(rd()); // Seed the generator
    std::shuffle(shuffled.begin(), shuffled.end(), gen);

    return std::vector<int>(shuffled.begin(), shuffled.begin() + t);

}

//threshold panta iso me 1 take into account unfiltered vectors which means they have all filters 
std::map<int, int> Vamana::Filtered_Find_Medoid(std::vector<Data_Point> dataset, std::unordered_set<int> filters, int threshold){

    std::map<int, int> M_map;
    std::map<int, int> T_map;
    
    int counter = 0;
    for(Data_Point point: dataset){
        T_map[counter] = 0;
        counter++;
    }

    for( const int &filter: filters ){

        std::vector<int> F_ids = get_ids(dataset, filter);
        std::vector<int> random_points = select_random_elements(F_ids, threshold);

        int min = std::numeric_limits<int>::max();
        int min_point = 0;

        for(int point: random_points){
            if(T_map[point] < min) {
                min = T_map[point];
                min_point = point;;
            }
        }

        M_map[filter] = min_point;
        T_map[min] = T_map[min] + 1; 
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
        minHeap.emplace(-1 * euclidean_distance(dataset[q].data_vector, dataset[node].data_vector), node);//* (-1) to make it min from max
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
            
            if(filterMin != -1){
                if(filterNode == -1){
                    if(filterQ != filterMin){
                        continue;
                    }
                }else if(filterQ == -1){
                    if(filterNode != filterMin){
                        continue;
                    }
                }else if(filterNode == filterQ){
                    if(filterNode != filterMin){
                        continue;
                    }
                }
            }
    
            double dis1 = euclidean_distance(dataset[minNode].data_vector, dataset[node].data_vector);
            double dis2 = euclidean_distance(dataset[q].data_vector, dataset[node].data_vector);

            if((a * dis1) <= dis2){
                V.erase(node);

                //remake heap for updated V
                minHeap = {};
                for(int i : V){
                    minHeap.emplace(-1 * euclidean_distance(dataset[q].data_vector, dataset[i].data_vector), i);
                }
            }
        }

    }

    //update new neighboors
    G.get_node(q)->neighbors = std::move(out);
}


RRGraph Vamana::Filtered_Vamana_Index(FilteredDataset dataset_obj, int L, int R, float a){

    std::unordered_set<int> filters_set = dataset_obj.get_filter_set();
    std::vector<Data_Point> dataset = dataset_obj.get_dataset();

    
    //create an empty graph
    RRGraph graph(R);
    graph.create_Rregular_empty_graph(dataset);

    //get the size of the dataset()
    int N = dataset.size();

    //get the random permutation as a starting 
    std::vector<int> perm = std::move(get_random_permutation(N)); 

    std::map<int, int> filter_map = Filtered_Find_Medoid(dataset, filters_set, 1);

    std::vector<int> visited;

    for(int i = 0; i<N; i++){

        cout << i << endl;

        // int starting_filter =  dataset[perm[i]].categorical;

        // int starting_node = filter_map[starting_filter];

        LVPair greedy_result = FilteredGreedySearch(graph, filter_map, perm[i], 0, L, filters_set, dataset);

        visited.insert(visited.end(), greedy_result.second.begin(), greedy_result.second.end());

        //run filtered robust prune
        FilteredRobustPruning(graph, perm[i], greedy_result.second, a, R, dataset_obj);

        //get the neighbors of perm[i]
        std::vector<int> perm_i_neighbours = graph.get_node(perm[i])->neighbors;

        for(int j: perm_i_neighbours){
            
            (graph.get_node(j)->neighbors).push_back(perm[i]);

            int n = graph.get_node(j)->neighbors.size();

            std::unordered_set<int> N_out_j(graph.get_node(j)->neighbors.begin(), graph.get_node(j)->neighbors.end());

            if(n > R){
                //run filtered robust prune
                FilteredRobustPruning(graph, j, N_out_j, a, R, dataset_obj);
            }

        }

    }

    for(int i = 0;i<graph.get_nodes_num();i++){
        if(graph.get_node(i)->neighbors.size() > R) graph.get_node(i)->neighbors.resize(R);
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
    //int medoid_index = 8736;

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

template <typename type>
GraphCollection Vamana::StitchedVamana(std::vector<std::vector<type>> dataset, std::unordered_set<int> filters, int Lsmall, int Rsmall, int Rstitched, float a){

    GraphCollection collectionOfGraphs;
    std::vector<std::vector<std::vector<type>>> PointsByFilter;
    int flag = 0;

    //create vector of points for each filter
    for(int f : filters.size()){
        for(int i = 0; i < dataset.size(); i++){
            if(dataset[i].filter == f){
                PointsByFilter[f].push_back(dataset[i]);
            }else if(dataset[i].filter == -1){
                PointsByFilter[f].push_back(dataset[i]);
                flag++;
            }
        }
    }

    //create subgraphs with vamana for each filter
    if(flag){
        int R = Rsmall;
    }esle{
        int R = Rstitched;
    }

    for(int f : filters.size()){
        RRGraph subGraph = Vamana::Vamana_Index(PointsByFilter[f], Lsmall, R, a);
        collectionOfGraphs.push_back(f, subGraph);
    } 


    if(flag){
        RRGraph stitchedGraph = RRGraph::create_Rregular_empty_graph(dataset);
        stitchedGraph = stitch_graph();
        collectionOfGraphs.clear();

        for (int q = 0; q < dataset.size(); q++) {
            FilteredRobustPruning(smallG, q, stitchedGraph.get_node(dataset[q])->neigbors, alpha, Rstitched, dataset);
        }
        collectionOfGraphs.push_back(0,stitchedGraph);
    }

    return collectionOfGraphs;

}

RRGraph RRGraph::stitch_graph(std::vector<std::vector<int, RRGraph>> G){

}

//function to get the recall 
double Vamana::Get_Recall(std::vector<int> vec1, std::vector<int> vec2){

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