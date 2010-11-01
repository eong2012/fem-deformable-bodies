#include "Face.h"

Face::Face(){

oppositeFaceInd = -1;
}
Face::~Face(){}

void Face::setNormal(Vector3<float> n)
{
    normal = n;
}
void Face::setEdgeInd(unsigned int e)
{
    edgeInd = e;

}
Vector3<float> Face::getNormal()
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