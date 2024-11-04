#ifndef GRAPH_H
#define GRAPH_H


#include <iostream>
#include <vector>


typedef struct{

    int node_id;
    std::vector<int> neighbors; 

}Node;

class RRGraph{

    private:
        //int a;
        int R;
        //int L;
        int nodes_num;
        std::vector<Node *> adj_list;

    public:
        RRGraph(int k);
        //RRGraph(const RRGraph& obj); //copy constructor
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
        
};


#endif //GRAPH_H
