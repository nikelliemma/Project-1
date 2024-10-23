#ifndef VAMANA_H
#define VAMANA_H

#include <iostream>
#include <vector>
#include <set>


#include "dataset.h"
#include "Graph.h"

typedef std::pair<std::vector<int>, std::set<int> > LVPair;

class Vamana{

    private:

        std::string filename;
        RRGraph * graph;
        
    public:

        Vamana();
        ~Vamana();
        template <typename Type>
        LVPair GreedySearch(int starting_node, int query, int k, int L, std::vector<vector<Type> > dataset);
        void RobustPruing();
        void VamanaIndex();

};



#endif