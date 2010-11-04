

#include "../Math/Vector3.h"

#ifndef FACE_H
#define FACE_H

class Face{
public:
    Face();
    ~Face();

    void setNormal(Vector3<float> n);
    void setEdgeInd(unsigned int e);
	void setOppositeFaceInd(int f);
    Vector3<float> getNormal();
    unsigned int getEdgeInd();
	int getOppositeFaceInd();


private:
    Vector3<float> normal;
    unsigned int edgeInd;
	int oppositeFaceInd;

};

#endif
