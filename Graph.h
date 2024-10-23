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
        int R;
        std::vector<Node *> adj_list;

    public:
        RRGraph(int k);
        ~RRGraph();
        int get_R();
        void set_R(int R);

        template <typename Type>
        void create_Rregular_graph(std::vector<std::vector<Type> > dataset);
        
};


#endif 