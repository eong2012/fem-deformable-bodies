#include "HalfEdge.h"
HalfEdge::HalfEdge(){}
HalfEdge::~HalfEdge(){}
unsigned int HalfEdge::getVertexInd()
{
    return vertexInd;
}
unsigned int HalfEdge::getFaceInd()
{
    return faceInd;
}
unsigned int HalfEdge::getPairInd()
{
    return pairInd;
}
unsigned int HalfEdge::getNextInd()
{
    return nextInd;
}
unsigned int HalfEdge::getPrevInd()
{
    return prevInd;
}

void HalfEdge::setVertexInd(unsigned int v)
{
    vertexInd = v;
}
void HalfEdge::setFaceInd(unsigned int f)
{
    faceInd = f;
}
void HalfEdge::setPairInd(unsigned int p)
{
    pairInd = p;
}
void HalfEdge::setNextInd(unsigned int n)
{
    nextInd = n;
}
void HalfEdge::setPrevInd(unsigned int p)
{
    prevInd = p;
}
