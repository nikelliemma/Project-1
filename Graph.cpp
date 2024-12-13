#include <iostream>
#include <vector>
#include <algorithm>
#include <random> 
#include <cstdint>

#include <fstream>
#include <stdexcept>

#include "Graph.h"
#include "filtered_dataset.h"


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

//templated function to create a random R-regular directed graph
template <typename Type>
void RRGraph::create_Rregular_graph(std::vector<std::vector<Type> > dataset){

    if(this->R <= 0) return;

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

            //ensure that the edge doesnt already exist in the neighbors vector
            while(random_number == i || std::find(new_node->neighbors.begin(), new_node->neighbors.end(), random_number) != new_node->neighbors.end()){
                random_number = dis(gen);
            }

            new_node->neighbors.push_back(random_number);
        }

        this->adj_list.push_back(new_node);
    }
    
    return;
}

// Function to create an empty R-regular graph
void RRGraph::create_Rregular_empty_graph(std::vector<Data_Point> dataset){
    int size = dataset.size(); // Number of data points in the dataset

    // Clear the adjacency list to ensure no leftover data from previous calls
    this->adj_list.clear();
    
    // Create a new node for each data point
    for(int i = 0; i < size; i++){
        Node* new_node = new Node;  // Dynamically allocate a new Node
        new_node->node_id = i;      // Assign a unique ID to the node
        this->adj_list.push_back(new_node);  // Add the node to the adjacency list
    }
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

// // Function to write the graph to a binary file
// void RRGraph::write_to_binary_file(const std::string& filename){
//     std::ofstream outfile(filename, std::ios::binary);
//     if (!outfile) {
//         throw std::runtime_error("Unable to open file for writing: " + filename);
//     }

//     // Write R and nodes_num
//     outfile.write(reinterpret_cast<const char*>(&this->R), sizeof(this->R));
//     outfile.write(reinterpret_cast<const char*>(&this->nodes_num), sizeof(this->nodes_num));

//     // Write each node and its neighbors
//     for (Node* node : this->adj_list) {
//         // Write node ID
//         outfile.write(reinterpret_cast<const char*>(&node->node_id), sizeof(node->node_id));

//         // Write the number of neighbors
//         int neighbors_count = node->neighbors.size();
//         outfile.write(reinterpret_cast<const char*>(&neighbors_count), sizeof(neighbors_count));

//         // Write the neighbors
//         for (int neighbor : node->neighbors) {
//             outfile.write(reinterpret_cast<const char*>(&neighbor), sizeof(neighbor));
//         }
//     }

//     outfile.close();
// }

void RRGraph::write_to_binary_file(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        throw std::runtime_error("Unable to open file for writing: " + filename);
    }

    // Write R and nodes_num
    outfile.write(reinterpret_cast<const char*>(&this->R), sizeof(this->R));
    outfile.write(reinterpret_cast<const char*>(&this->nodes_num), sizeof(this->nodes_num));

    // Write each node and its neighbors
    for (Node* node : this->adj_list) {
        // Write node ID
        outfile.write(reinterpret_cast<const char*>(&node->node_id), sizeof(node->node_id));

        // Write the number of neighbors
        int neighbors_count = static_cast<int>(node->neighbors.size());
        outfile.write(reinterpret_cast<const char*>(&neighbors_count), sizeof(neighbors_count));

        // Write the neighbors
        if (neighbors_count > 0) {
            outfile.write(reinterpret_cast<const char*>(node->neighbors.data()),
                          neighbors_count * sizeof(int));
        }
    }

    outfile.close();
}




void RRGraph::read_from_binary_file(const std::string& filename) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Unable to open file for reading: " + filename);
    }

    // Clear the existing adjacency list
    for (Node* node : this->adj_list) {
        delete node;
    }
    this->adj_list.clear();

    // Read R and nodes_num
    infile.read(reinterpret_cast<char*>(&this->R), sizeof(this->R));
    if (!infile) {
        throw std::runtime_error("Error reading R from file");
    }
    infile.read(reinterpret_cast<char*>(&this->nodes_num), sizeof(this->nodes_num));
    if (!infile) {
        throw std::runtime_error("Error reading nodes_num from file");
    }

    //std::cout << "Read R = " << this->R << ", nodes_num = " << this->nodes_num << std::endl;

    // Read each node and its neighbors
    for (int i = 0; i < this->nodes_num; i++) {
        Node* new_node = new Node;

        // Read node ID
        infile.read(reinterpret_cast<char*>(&new_node->node_id), sizeof(new_node->node_id));
        if (!infile) {
            throw std::runtime_error("Error reading node ID from file");
        }

        // Read the number of neighbors
        int neighbors_count = 0;
        infile.read(reinterpret_cast<char*>(&neighbors_count), sizeof(neighbors_count));
        if (!infile) {
            throw std::runtime_error("Error reading neighbors_count from file");
        }

        //std::cout << "Reading Node " << new_node->node_id << " with " << neighbors_count << " neighbors" << std::endl;

        // Read the neighbors
        if (neighbors_count > 0) {
            new_node->neighbors.resize(neighbors_count);
            infile.read(reinterpret_cast<char*>(new_node->neighbors.data()), neighbors_count * sizeof(int));
            if (!infile) {
                throw std::runtime_error("Error reading neighbors from file");
            }
        }

        this->adj_list.push_back(new_node);
    }

    infile.close();

}



// void RRGraph::read_from_binary_file(const std::string& filename) {
//     std::ifstream infile(filename, std::ios::binary);
//     if (!infile) {
//         throw std::runtime_error("Unable to open file for reading: " + filename);
//     }
    
//     // Clear the existing adjacency list
//     for (Node* node : this->adj_list) {
//         delete node;
//     }
//     this->adj_list.clear();

//     // Read R and nodes_num
//     infile.read(reinterpret_cast<char*>(&this->R), sizeof(this->R));
//     infile.read(reinterpret_cast<char*>(&this->nodes_num), sizeof(this->nodes_num));

//     // Read each node and its neighbors
//     for (int i = 0; i < this->nodes_num; i++) {
//         Node* new_node = new Node;

//         // Read node ID
//         infile.read(reinterpret_cast<char*>(&new_node->node_id), sizeof(new_node->node_id));

//         // Read the number of neighbors
//         int neighbors_count = 0;
//         infile.read(reinterpret_cast<char*>(&neighbors_count), sizeof(neighbors_count));

//         // Read the neighbors
//         if (neighbors_count > 0) {
//             new_node->neighbors.resize(neighbors_count);
//             infile.read(reinterpret_cast<char*>(new_node->neighbors.data()),
//                         neighbors_count * sizeof(int));
//         }

//         this->adj_list.push_back(new_node);
//     }

//     infile.close();
// }


// // Function to read the graph from a binary file
// void RRGraph::read_from_binary_file(const std::string& filename) {
//     std::ifstream infile(filename, std::ios::binary);
//     if (!infile) {
//         throw std::runtime_error("Unable to open file for reading: " + filename);
//     }

//     // Clear the existing adjacency list
//     for (Node* node : this->adj_list) {
//         delete node;
//     }
//     this->adj_list.clear();

//     // Read R and nodes_num
//     infile.read(reinterpret_cast<char*>(&this->R), sizeof(this->R));
//     infile.read(reinterpret_cast<char*>(&this->nodes_num), sizeof(this->nodes_num));

//     // Read each node and its neighbors
//     for (int i = 0; i < this->nodes_num; i++) {
//         Node* new_node = new Node;

//         // Read node ID
//         infile.read(reinterpret_cast<char*>(&new_node->node_id), sizeof(new_node->node_id));

//         // Read the number of neighbors
//         int neighbors_count;
//         infile.read(reinterpret_cast<char*>(&neighbors_count), sizeof(neighbors_count));

//         // Read the neighbors
//         new_node->neighbors.resize(neighbors_count);
//         for (int j = 0; j < neighbors_count; j++) {
//             infile.read(reinterpret_cast<char*>(&new_node->neighbors[j]), sizeof(new_node->neighbors[j]));
//         }

//         this->adj_list.push_back(new_node);
//     }

//     infile.close();
// }




template void RRGraph::create_Rregular_graph<float>(std::vector<std::vector<float> > dataset);
template void RRGraph::create_Rregular_graph<int>(std::vector<std::vector<int> > dataset);
template void RRGraph::create_Rregular_graph<unsigned char>(std::vector<std::vector<unsigned char> > dataset);
