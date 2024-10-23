#include <iostream>
#include "dataset.h"

using namespace std;


int main(void){


    Dataset<float> d;

    d.set_filename(".fvecs"); //add the fvecs file

    d.set_type(FLOAT); //set type FLOAT/INTEGER/UN_CHAR

    d.read_dataset();

    int dim = d.get_dim();

    //print only the vector 0 for confirmation
    for(int i = 0; i<d.get_dim(); i++){
        cout << d.get_data(0, i) << " ";
    }
    cout << endl;

    cout << "vector dimension : " << d.get_dim() << endl;

    cout << "total vector number : " << d.get_vectors_num() << endl;


    return 0;
}