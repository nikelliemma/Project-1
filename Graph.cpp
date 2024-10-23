#include <iostream>
#include <vector>
#include <algorithm>
#include <random> 

#include "Graph.h"


RRGraph::RRGraph(int k){
    //this->groundtruth_flag = false;
    this->R = k; 
    return;
}

RRGraph::~RRGraph(){

    for(Node* node : this->adj_list){
        delete node;
    }

    this->adj_list.clear();
}

int RRGraph::get_R() {
    return this->R;
}

void RRGraph::set_R(int R) {
    this->R = R;
}



//NOT FINISHED
template <typename Type>
void RRGraph::create_Rregular_graph(std::vector<std::vector<Type> > dataset){

    int vecs_num = dataset.size();
    
    for(int i = 0; i<vecs_num; i++){

        Node * new_node = new Node;
        new_node->node_id = i;

        std::random_device rd;  
        std::mt19937 gen(rd()); 
        std::uniform_int_distribution<> dis(0, vecs_num - 1);

        for(int j = 0; j < R; j++){

            int random_number = dis(gen);

            while(random_number == i || std::find(new_node->neighbors.begin(), new_node->neighbors.end(), random_number) != new_node->neighbors.end()){
                random_number = dis(gen);
            }

            new_node->neighbors.push_back(random_number);
        }

        this->adj_list.push_back(new_node);
    }
    
    return;
}
