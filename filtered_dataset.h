#ifndef FILTERED_DATASET_H
#define FILTERED_DATASET_H

#include <iostream>
#include <vector>
#include <unordered_set>


typedef struct{
    std::vector<float> data_vector;          
    int categorical;
    float timestamp;
}Data_Point;

//class that stores the dataset
//inlcudes a method that reads the dataset from the binary file
class FilteredDataset{

    private:
    //private members of the Dataset class 
        std::string file_path; //contains the filepath
        std::vector<Data_Point> dataset;
        int vectors_num;
        int dimension;
        std::unordered_set<int> filter_set;


    public:
    //class method definitions
        FilteredDataset(); //constructor
        ~FilteredDataset(); //destructor

        //getters - setters
        void set_vectors_num(int vectors_num);
        int get_vectors_num();
        void set_dimension(int dimension);
        int get_dimension();
        std::string get_filepath();
        void set_filepath(std::string filepath);
        Data_Point get_data_point(int index);
        std::vector<Data_Point> get_dataset();
        std::unordered_set<int> get_filter_set();

        //read dataset and query set methods
        void read_Dataset();
        void read_Query_set();

};

#endif
