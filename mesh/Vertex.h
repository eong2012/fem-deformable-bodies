
#include "../Math/Vector3.h"

#ifndef VERTEX_H
#define VERTEX_H

class Vertex{
public:
    Vertex();
    ~Vertex();

    void setPosition(Vector3<float> p);
    void setEdgeInd(unsigned int e);
    Vector3<float> getPosition();
    unsigned int getEdgeInd();


private:
    Vector3<float> position;
    unsigned int edgeInd;

};
#endif
