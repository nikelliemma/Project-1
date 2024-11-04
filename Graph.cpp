#include <iostream>
#include <vector>
#include <algorithm>
#include <random> 

#include "Graph.h"


RRGraph::RRGraph(int R){
    this->R = R; 
    this->nodes_num = 0;
    return;
}

RRGraph::~RRGraph(){

}

int RRGraph::get_R(){
    return this->R;
}

void RRGraph::set_R(int R){
    this->R = R;
}

Node * RRGraph::get_node(int index){
    return this->adj_list[index];
}

int RRGraph::get_nodes_num(){
    return this->nodes_num;
}

void RRGraph::set_nodes_num(int num){
    this->nodes_num = num;
}


std::vector<Node *> RRGraph::get_graph(){
    return this->adj_list;
}

//templated function to create an R-regular directed graph
template <typename Type>
void RRGraph::create_Rregular_graph(std::vector<std::vector<Type> > dataset){

    int vecs_num = dataset.size();
    set_nodes_num(vecs_num);
    
    for(int i = 0; i<vecs_num; i++){

        Node * new_node = new Node;
        new_node->node_id = i;

        std::random_device rd;  
        std::mt19937 gen(rd()); 
        std::uniform_int_distribution<> dis(0, vecs_num - 1);

        for(int j = 0; j < this->R; j++){

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

//utlity function to print the graph
void RRGraph::print_graph(){

    for (Node* node : adj_list){

        std::cout << "Node " << node->node_id << ": ";

        for(int neighbor : node->neighbors){
            std::cout << neighbor << " ";
        }

        std::cout << std::endl;
    }
}



template void RRGraph::create_Rregular_graph<float>(std::vector<std::vector<float> > dataset);
template void RRGraph::create_Rregular_graph<int>(std::vector<std::vector<int> > dataset);
template void RRGraph::create_Rregular_graph<unsigned char>(std::vector<std::vector<unsigned char> > dataset);