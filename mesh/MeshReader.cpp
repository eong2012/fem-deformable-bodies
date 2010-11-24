#include "MeshReader.h"

MeshReader::MeshReader(){}
MeshReader::~MeshReader(){}

vector<unsigned int> MeshReader::readTetras(char *filename){
    vector<unsigned int> indices;

    ifstream infile;
    infile.open(filename);

    string trash;
    getline(infile, trash);
    int tmp, ind0, ind1, ind2, ind3;

    while (infile.good()){
        //infile >> tmp;
        infile >> ind0;
        infile >> ind1;
        infile >> ind2;
        infile >> ind3;
        indices.push_back(ind0);
        indices.push_back(ind1);
        indices.push_back(ind2);
        indices.push_back(ind3);
        getline(infile, trash);

    }
    cout << indices.size() << endl;
    infile.close();
    return indices;

    }

vector<arma::Mat<double> > MeshReader::readVertices(char *filename){
    vector<arma::Mat<double> > vertices;

    ifstream infile;
    infile.open(filename);

    string trash;
    getline(infile, trash);

    double tmp, x, y, z;
    arma::Mat<double> tmpVert;

    while (infile.good()){
        //infile >> tmp;
        infile >> x;
        infile >> y;
        infile >> z;
        tmpVert << x << y << z;
        vertices.push_back(tmpVert);
        getline(infile, trash);
    }


    infile.close();
    return vertices;

    }
