#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>

#include "filtered_dataset.h"


//--------------Constructor/Dstructor--------------//

FilteredDataset::FilteredDataset(){
    this->file_path = "";
    this->dimension = 100;
    this->vectors_num = 0;
}

FilteredDataset::~FilteredDataset(){

}


//-------------------Get Methods-------------------//

int FilteredDataset::get_vectors_num(){
    return this->vectors_num;
}

int FilteredDataset::get_dimension(){
    return this->dimension;
}

std::string FilteredDataset::get_filepath(){
    return this->file_path;
}

std::vector<Data_Point> FilteredDataset::get_dataset(){
    return this->dataset;
}

Data_Point FilteredDataset::get_data_point(int index){

    return this->dataset[index];
}

std::unordered_set<int> FilteredDataset::get_filter_set(){
    return this->filter_set;
}

int FilteredDataset::get_query_type(int index){
    return this->dataset[index].timestamp;
}

//-------------------Set Methods-------------------//

void FilteredDataset::set_dimension(int dimension){
    this->dimension = dimension;
}

void FilteredDataset::set_vectors_num(int vectors_num){
    this->vectors_num = vectors_num;
}

void FilteredDataset::set_filepath(std::string filepath){
    this->file_path = filepath;
}



//---------------Read Dataset Method---------------//

void FilteredDataset::read_Dataset(){

    std::ifstream file(this->get_filepath(), std::ios::binary);
    if(!file.is_open()){
        throw std::runtime_error("Failed to open the binary file.");
    }

    uint32_t vectors_num;
    file.read(reinterpret_cast<char *>(&vectors_num), sizeof(uint32_t));
    this->set_vectors_num(vectors_num);

    for(int i=0; i<vectors_num; i++){
        Data_Point point;

        float categorical;
        file.read(reinterpret_cast<char *>(&categorical), sizeof(float));
        point.categorical = static_cast<int>(categorical);
        this->filter_set.emplace(categorical);

        float timestamp;
        file.read(reinterpret_cast<char *>(&timestamp), sizeof(float));
        point.timestamp = timestamp;

        point.data_vector.resize(this->get_dimension());
        file.read(reinterpret_cast<char *>(point.data_vector.data()), sizeof(float) * this->get_dimension());

        dataset.push_back(point);
    }

    return;
}

void FilteredDataset::read_Query_set(){

    // Open the binary file for reading
    std::ifstream file(this->get_filepath(), std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the binary file.");
    }

    // Read the number of queries from the file
    uint32_t num_queries;
    file.read(reinterpret_cast<char*>(&num_queries), sizeof(uint32_t));

    // Read each query from the file
    for (uint32_t i = 0; i < num_queries; i++) {
        Data_Point query;

        // Read the first four dimensions
        float query_type;
        file.read(reinterpret_cast<char*>(&query_type), sizeof(float));
        query.timestamp = static_cast<int>(query_type);

        float categorical_value;
        file.read(reinterpret_cast<char*>(&categorical_value), sizeof(float));
        query.categorical = static_cast<int>(categorical_value);

        float timestamp_lower;
        file.read(reinterpret_cast<char*>(&timestamp_lower), sizeof(float));

        float timestamp_upper;
        file.read(reinterpret_cast<char*>(&timestamp_upper), sizeof(float));

        // Read the rest of the 100 dimensions (query vector)
        query.data_vector.resize(this->get_dimension());
        file.read(reinterpret_cast<char*>(query.data_vector.data()), sizeof(float) * this->get_dimension());

        // Add the query to the query set
        //if(query_type != 0 && query_type != 1) continue;
        dataset.push_back(query);
        
    }

    return;
}


