#ifndef DATASET_H
#define DATASET_H

#include <iostream>
#include <vector>

using namespace std; 


// #define cast(ptr, type) (*static_cast<type*>(ptr))


const string BVECS_FORMAT = "bvecs";
const string IVECS_FORMAT = "ivecs";
const string FVECS_FORMAT = "fvecs";

const string INTEGER = "integer";
const string FLOAT = "float";
const string UN_INTEGER = "un_integer";

template <typename Type>

class Dataset{

    private:   
        string type; 
        string filename;
        int dim;
        vector<vector<Type> > dataset;
        int vectors_num;

    public:

        Dataset();
        ~Dataset();
        string get_filename();
        void set_filename(string filename);
        int get_dim();
        void set_dim(int dim);
        int get_vectors_num();
        void set_vectors_num(int vec_num);
        string extract_format(string &file);
        void read_dataset();
        Type get_data(int vec_num, int index);
        void set_type(string type);
        string get_type();
        vector<Type> get_vector(int index);

        template <typename TempType>
        void set_dataset(std::vector<std::vector<TempType> >& source);
        vector<vector<Type> > get_dataset();
        
};


#endif //DATASET_H
