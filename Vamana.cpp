#include <iostream>
#include <vector>
#include <set>
#include <cmath>



#include "Vamana.h"
#include "Graph.h"


//constructor
Vamana::Vamana(){

}

//destructor
Vamana::~Vamana(){

}

bool check_loop_cond(std::vector<int> L, std::set<int> V){

    for(int i = 0; i<L.size(); i++){
        if(V.find(L[i]) != V.end()) return true;
    }

    return false;
}


//function to compute Euclidean distance between two vectors
template <typename type>
double euclidean_distance(const std::vector<type>& vec1, const std::vector<type>& vec2){

    double dist = 0.0;
    for(int i = 0; i < vec1.size(); ++i){

        double diff = vec1[i] - vec2[i];
        dist += diff * diff;

    }

    return std::sqrt(dist);
}


template <typename type>
int get_medoid(std::vector<std::vector<type> > datastet){
    
}


//NOT FINISHED!
template <typename type>
LVPair Vamana::GreedySearch(int starting_node, int query, int k, int L, std::vector<std::vector<type> > datastet){

    //result vector contains the NNs of the query
    //visited set contains the visited nodes 
    std::vector<int> result;
    std::set<int> visited;


    result.push_back(starting_node);

    //start the greedy search loop
    while(check_loop_cond(result, visited)){


        //compute the euclidian distance
        

        //retain the L closest points to the query
        std::sort(result.begin(), result.end());
        if(result.size() > L){
            result.resize(L);
        }
    }


    //get the k nearest neighbors
    std::sort(result.begin(), result.end());
    result.resize(k);

    //return the NNs and the visited nodes set
    return {result, visited};
}




void Vamana::RobustPruing(){


//-------------------------------


// EMMA


//-------------------------------

    return;
}


void Vamana::VamanaIndex(){



    return;
}