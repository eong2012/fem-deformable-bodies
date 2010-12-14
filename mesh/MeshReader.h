#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Vertex.h"
#include "armadillo"


using namespace std;

class MeshReader{
    public:
    MeshReader();
    ~MeshReader();

    vector<arma::Mat<double> > readVertices(char *filename);
    vector<unsigned int> readTetras(char *filename);


    private:




    };
