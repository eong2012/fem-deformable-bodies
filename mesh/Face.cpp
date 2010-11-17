#include "Face.h"

Face::Face(){

oppositeFaceInd = -1;
}
Face::~Face(){}

void Face::setNormal(arma::Mat<double> n)
{
    normal = n;
}
void Face::setEdgeInd(unsigned int e)
{
    edgeInd = e;

}
arma::Mat<double> Face::getNormal()
{
    return normal;
}
unsigned int Face::getEdgeInd()
{
    return edgeInd;
}

void Face::setOppositeFaceInd(int f) {

	this->oppositeFaceInd = f;
}

int Face::getOppositeFaceInd() {

return this->oppositeFaceInd;
}
