#ifndef VAMANA_H
#define VAMANA_H

#include <iostream>
#include <vector>
#include <unordered_set>


#include "dataset.h"
#include "Graph.h"

typedef std::pair<std::vector<int>, std::unordered_set<int> > LVPair;

typedef std::vector<std::vector<int, RRGraph>> GraphCollection; //new, might change

class Vamana{

    private:

        int L;
        int R;
        int alpha;
        std::string filename;
        RRGraph vamana_index;
        
    public:

        Vamana(int R, int L, int alpha);
        ~Vamana();
        int get_L();
        int get_R();
        int get_alpha();
        void set_L(int L);
        void set_R(int R);
        void set_alpha(int alpha);
        void set_index(RRGraph graph);
        RRGraph get_index();
        template <typename type>
        double euclidean_distance(const std::vector<type>& vec1, const std::vector<type>& vec2);
        template <typename Type>
        LVPair GreedySearch(RRGraph graph, int starting_node, int query, int k, int L, std::vector<vector<Type> > dataset);
        template <typename Type>
        LVPair GreedySearch(RRGraph graph, int starting_node, std::vector<Type> q_vec, int k, int L, std::vector<vector<Type> > dataset);
        template <typename type>
        RRGraph Vamana_Index(std::vector<std::vector<type> > dataset, int L, int R, float a);
        template <typename type>
        RRGraph StitchedVamana(std::vector<std::vector<type> > dataset, std::unordered_set<int> filters, int Lsmall, int Rsmall, int Rstitched, float a);
        template <typename type>
        void RobustPruning(RRGraph graph, int query, std::unordered_set<int> V_set, float a, int R, std::vector<std::vector<type> > dataset);
        template <typename type>
        void FilteredRobustPruning(RRGraph G, int q, std::unordered_set<int> V, float a, int R, std::vector<std::vector<type>> dataset);
        template <typename type>
        int find_medoid(std::vector<std::vector<type> > dataset);
        double Get_Recall(std::vector<int> vec1, std::vector<int> vec2);
        void create_vamana_index(std::string filepath, int L, int R, int alpha);

};


#endif