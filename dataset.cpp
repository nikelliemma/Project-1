#include <iostream>
#include <vector>
#include <fstream>

#include "dataset.h"


using namespace std;



template <typename Type>
Dataset<Type>::Dataset(){
    this->type = "";
    this->filename = "";
    this->dim = 0;
    this->vectors_num = 0;
}

template <typename Type>
Dataset<Type>::~Dataset(){ }

template <typename Type>
string Dataset<Type>::get_filename(){ return this->filename; }

template <typename Type>
void Dataset<Type>::set_filename(string filename){ this->filename = filename; }

template <typename Type>
int Dataset<Type>::get_dim(){ return this->dim; }

template <typename Type>
void Dataset<Type>::set_dim(int dim){ this->dim = dim; }

template <typename Type>
void Dataset<Type>::set_vectors_num(int vec_num){ this->vectors_num = vec_num; }

template <typename Type>
int Dataset<Type>::get_vectors_num(){ return this->vectors_num; }

template <typename Type>
string Dataset<Type>::get_type(){ return this->type; }

template <typename Type>
void Dataset<Type>::set_type(string type){ this->type = type; }



template <typename Type>
string Dataset<Type>::extract_format(string &file){

    string format;
    int len = file.length();

    int pos = file.find_last_of("/");
    string file_string = file.substr(pos + 1);

    int index = 0;
    while(file[index] != '.'){ index++; }

    if(index == len - 1) return "";

    format = file.substr(index + 1);

    return format;

}


int read_dimension(std::ifstream &file){

    int dimension;

    file.read(reinterpret_cast<char *>(&dimension), sizeof(dimension));
    
    // if(file.fail()){
    //     cerr << "Error while reading the dimension of the vectors." << endl;
    //     return -1;
    // }

    //this->dim = dimension;

    return dimension;
}

vector<vector<float> > read_fvecs(string filename){

    vector<vector <float> > dataset;

    std::ifstream file;
    file.open(filename, std::ios::binary);
    
    if(!file.is_open()){
        cerr << "Error while opening file: " << filename << endl;
        return dataset;
    }

    int dim;

    while(!file.eof()){

        dim = read_dimension(file);

        if(!file) break;
        
        if(dim <= 0){
            cout << "Error: Zero or Negative dimension!" << endl;
            file.close();
            return dataset; //return an empty dataset
        }

        std::size_t size_to_read = dim * sizeof(float);

        vector<float> util_vec(dim);

        file.read(reinterpret_cast<char *>(util_vec.data()), size_to_read);

        dataset.push_back(util_vec);

    }

    file.close();

    return dataset;

}

vector<vector<int> > read_ivecs(string filename){

    vector<vector <int> > dataset;

    std::ifstream file;
    file.open(filename, std::ios::binary);
    
    if(!file.is_open()){
        cerr << "Error while opening file: " << filename << endl;
        file.close();
        return dataset;
    }

    int dim;

    while(!file.eof()){

        dim = read_dimension(file);

        if(!file) break;

        if(dim <= 0){
            cout << "Error: Zero or Negative dimension!" << endl;
            return dataset; //return an empty dataset
        }

        std::size_t size_to_read = dim * sizeof(float);

        vector<int> util_vec(dim);

        file.read(reinterpret_cast<char *>(util_vec.data()), size_to_read);

        dataset.push_back(util_vec);

    }

    file.close();

    return dataset;

}

vector<vector<unsigned char> > read_bvecs(string filename){

    vector<vector <unsigned char> > dataset;

    std::ifstream file;
    file.open(filename, std::ios::binary);
    
    if(!file.is_open()){
        cerr << "Error while opening file: " << filename << endl;
        return dataset;
    }

    int dim;

    while(!file.eof()){

        dim = read_dimension(file);

        if(!file) break;

        if(dim <= 0){
            cout << "Error: Zero or Negative dimension!" << endl;
            file.close();
            return dataset; //return an empty dataset
        }

        std::size_t size_to_read = dim * sizeof(float);

        vector<unsigned char> util_vec(dim);

        file.read(reinterpret_cast<char *>(util_vec.data()), size_to_read);

        dataset.push_back(util_vec);

    }

    file.close();

    return dataset;

}


template <typename Type>
template <typename TempType>
void Dataset<Type>::set_dataset(std::vector<std::vector<TempType> >& source){

    this->dataset.clear();  //clear any existing data
    this->dataset.reserve(source.size());  //reserve space for efficiency

    for(auto& inner_vec : source){

        std::vector<Type> temp_vec;
        temp_vec.reserve(inner_vec.size());  //reserve space for inner vectors

        for(TempType& element : inner_vec){
            temp_vec.push_back(static_cast<Type>(element));   //cast
        }

        this->dataset.push_back(std::move(temp_vec));  //move temp_vec to data
    }   

    return;
}


template <typename Type>
void Dataset<Type>::read_dataset(){

    string file_str = this->filename;
    string format = extract_format(file_str);

    if(format == BVECS_FORMAT){
        //read .bvecs file format
        vector<vector<unsigned char> > temp_dataset = read_bvecs(file_str);
        set_dataset(temp_dataset);
        this->type = UN_CHAR;
        this->dim = temp_dataset[0].size();
        this->vectors_num = temp_dataset.size();
    }
    else if(format == IVECS_FORMAT){
        //read .ivecs file format
        vector<vector<int> > temp_dataset = read_ivecs(file_str);
        set_dataset(temp_dataset);
        this->type = INTEGER;
        this->dim = temp_dataset[0].size();
        this->vectors_num = temp_dataset.size();
    }
    else if(format == FVECS_FORMAT){
        //read .ivecs file format
        vector<vector<float> > temp_dataset = read_fvecs(file_str);
        set_dataset(temp_dataset);
        this->type = FLOAT;
        this->dim = temp_dataset[0].size();
        this->vectors_num = temp_dataset.size();
    }
    else{
        throw std::runtime_error("Unsupported format");
    }

    return;
}

template <typename Type>
Type Dataset<Type>::get_data(int vec_num, int index){
    return this->dataset[vec_num][index];
}

template <typename Type>
vector<Type> Dataset<Type>::get_vector(int index){
    return this->dataset[index];
}

//explicit template instantiation
template class Dataset<float>;
template class Dataset<int>;
template class Dataset<unsigned char>;
