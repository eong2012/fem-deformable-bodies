
#include "armadillo"

#ifndef VERTEX_H
#define VERTEX_H

class Vertex{
public:
    Vertex();
    ~Vertex();

    void setPosition(arma::Mat<double> p);
    void setEdgeInd(unsigned int e);
    arma::Mat<double> getPosition();
    unsigned int getEdgeInd();


private:
    arma::Mat<double> position;
    unsigned int edgeInd;

};
#endif
