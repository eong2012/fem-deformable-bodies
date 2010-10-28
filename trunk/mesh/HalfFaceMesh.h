#include "HalfFace.h"
#include "Face.h"
#include "Vertex.h"
#include "HalfEdge.h"
#include <vector>
using namespace std;

class HalfFaceMesh{
public:
    HalfFaceMesh();
    ~HalfFaceMesh();

    void AddTetrahedron();
    void AddHalfFace(vector<Vector3<float> > vertices);
    void AddHalfEdgePair(unsigned int vertexIndex1,unsigned int vertexIndex2,unsigned int &edgeIndex1,unsigned int &edgeIndex2);
    bool AddVertex(Vector3<float> vertexPos, unsigned int &index);

    void setTetraGeometry();

private:

    vector<Vertex> mVertices;
    vector<HalfEdge> mHalfEdges;
    vector<Face> mFaces;
    vector<HalfFace> mHalfFaces;



};
