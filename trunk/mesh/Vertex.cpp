#include "Vertex.h"
Vertex::Vertex(){}
Vertex::~Vertex(){}

void Vertex::setPosition(arma::Mat<double> p)
{
   position = p;
}
void Vertex::setEdgeInd(unsigned int e)
{
    edgeInd = e;
}
arma::Mat<double> Vertex::getPosition()
{
    return position;
}
unsigned int Vertex::getEdgeInd()
{
    return edgeInd;
}
