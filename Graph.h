#ifndef GRAPH_H
#define GRAPH_H


#include <iostream>
#include <vector>

#include "filtered_dataset.h"

//this class represents an R-regular graph in adjacency list representation,
//where the adjacency list is a vector of a struct Node * object, this Node struct object
//contains the node ID and a vectors of the neighbors/outgoing edges of that node. 


//node struct, it contains the node id and the neighbors list (the outgoing edges of the node)
typedef struct{

    int node_id;
    std::vector<int> neighbors; 

}Node;

class RRGraph{

    private:
        int R;
        int nodes_num;
        std::vector<Node *> adj_list;

    public:
        RRGraph(int R);
        ~RRGraph();
        int get_R();
        void set_R(int R);
        Node * get_node(int index);
        void print_graph();
        int get_nodes_num();
        void set_nodes_num(int num);
        std::vector<Node *> get_graph();

        template <typename Type>
        void create_Rregular_graph(std::vector<std::vector<Type> > dataset);
        void create_Rregular_empty_graph(std::vector<Data_Point> dataset);

        void write_to_binary_file(const std::string& filename);
        void read_from_binary_file(const std::string& filename);
        
};


#endif //GRAPH_H
