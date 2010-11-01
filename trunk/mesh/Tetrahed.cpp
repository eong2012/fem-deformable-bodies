#include "Tetrahed.h"
Tetrahed::Tetrahed(){}
Tetrahed::~Tetrahed(){}

vector<unsigned int> Tetrahed::getFaceInd()
{
	vector<unsigned int> temp;
	temp.push_back(this->faceInd1);
	temp.push_back(this->faceInd2);
	temp.push_back(this->faceInd3);
	temp.push_back(this->faceInd4);

    return temp;
}


void Tetrahed::setFaceInd1(unsigned int f)
{
    faceInd1 = f;
}

void Tetrahed::setFaceInd2(unsigned int f)
{
    faceInd2 = f;
}

void Tetrahed::setFaceInd3(unsigned int f)
{
    faceInd3 = f;
}

void Tetrahed::setFaceInd4(unsigned int f)
{
    faceInd4 = f;
}