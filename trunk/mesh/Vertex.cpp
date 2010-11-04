#include "Vertex.h"
Vertex::Vertex(){}
Vertex::~Vertex(){}

void Vertex::setPosition(Vector3<float> p)
{
   position = p;
}
void Vertex::setEdgeInd(unsigned int e)
{
    edgeInd = e;
}
Vector3<float> Vertex::getPosition()
{
    return position;
}
unsigned int Vertex::getEdgeInd()
{
    return edgeInd;
}
