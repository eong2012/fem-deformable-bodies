#include "HalfFace.h"
HalfFace::HalfFace(){}
HalfFace::~HalfFace(){}
unsigned int HalfFace::getVertexInd()
{
    return vertexInd;
}
unsigned int HalfFace::getEdgeInd()
{
    return edgeInd;
}
unsigned int HalfFace::getFaceInd()
{
    return faceInd;
}
unsigned int HalfFace::getMateInd()
{
    return mateInd;
}
unsigned int HalfFace::getNextInd()
{
    return nextInd;
}
unsigned int HalfFace::getPrevInd()
{
    return prevInd;
}
unsigned int HalfFace::getRadialInd()
{
    return radialInd;
}

void HalfFace::setVertexInd(unsigned int v)
{
    vertexInd = v;
}
void HalfFace::setEdgeInd(unsigned int e)
{
    edgeInd = e;
}
void HalfFace::setFaceInd(unsigned int f)
{
    faceInd = f;
}
void HalfFace::setMateInd(unsigned int m)
{
    mateInd = m;
}
void HalfFace::setNextInd(unsigned int n)
{
    nextInd = n;
}
void HalfFace::setPrevInd(unsigned int p)
{
    prevInd = p;
}
void HalfFace::setRadialInd(unsigned int r)
{
    radialInd = r;
}
