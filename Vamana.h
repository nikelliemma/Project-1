#ifndef VAMANA_H
#define VAMANA_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <map>

#include "dataset.h"
#include "filtered_dataset.h"
#include "Graph.h"


#define MAX_DIS (100000)

typedef std::pair<std::vector<int>, std::unordered_set<int> > LVPair;

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
        double euclidean_distance(const std::vector<int>& vec1, const std::vector<int>& vec2);
        template <typename type>
        double euclidean_distance(const std::vector<type>& vec1, const std::vector<type>& vec2, double min_dis);
        template <typename Type>
        LVPair GreedySearch(RRGraph graph, int starting_node, int query, int k, int L, std::vector<vector<Type> > dataset);
        template <typename Type>
        LVPair GreedySearch(RRGraph graph, int starting_node, std::vector<Type> q_vec, int k, int L, std::vector<vector<Type> > dataset);
        template <typename type>
        RRGraph Vamana_Index(std::vector<std::vector<type> > dataset, int L, int R, float a);
        template <typename type>
        void RobustPruning(RRGraph graph, int query, std::unordered_set<int> V_set, float a, int R, std::vector<std::vector<type> > dataset);
        template <typename type>
        int find_medoid(std::vector<std::vector<type> > dataset);
        int find_medoid_filtered(std::vector<Data_Point> dataset, int filter);
        double Get_Recall(std::vector<int> vec1, std::vector<int> vec2);
        double Get_Recall_Filtered(std::vector<int> knns, std::vector<int> vec2, FilteredDataset dataset, int k, int filter);
        void create_vamana_index(std::string filepath, int L, int R, int alpha);
        void FilteredRobustPruning(RRGraph G, int q, std::unordered_set<int> V, float a, int R, FilteredDataset f_dataset);
        LVPair FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset, std::vector<Data_Point> Q_dataset);
        LVPair FilteredGreedySearch(RRGraph graph, std::map<int, int> S_nodes, int query_vec, int k, int L, std::unordered_set<int> filters, std::vector<Data_Point> dataset);
        std::map<int, int> Filtered_Find_Medoid(std::vector<Data_Point> dataset, std::unordered_set<int> filters,int threshold);
        RRGraph Filtered_Vamana_Index(FilteredDataset dataset_obj, int L, int R, float a);
        RRGraph Filtered_Vamana_Index_Parallel(FilteredDataset dataset_obj, int L, int R, float a); 
        GraphCollection StitchedVamanaParallel(FilteredDataset dataset_obj, int Lsmall, int Rsmall, int Rstitched, float a);
        GraphCollection StitchedVamana(FilteredDataset, int , int , int , float);
        int find_medoid_f(std::vector<Data_Point> dataset);

};


#endif
