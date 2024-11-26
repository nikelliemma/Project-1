#include <iostream>
#include <vector>
#include <fstream>

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

    //if(index > this->get_vectors_num()) return ;

    return this->dataset[index];
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
    return;
}
