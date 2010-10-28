

#include "../Math/Vector3.h"
class Face{
public:
    Face();
    ~Face();

    void setNormal(Vector3<float> n);
    void setEdgeInd(unsigned int e);
    Vector3<float> getNormal();
    unsigned int getEdgeInd();


private:
    Vector3<float> normal;
    unsigned int edgeInd;

};

