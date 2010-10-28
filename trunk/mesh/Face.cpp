#include "Face.h"

Face::Face(){}
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
