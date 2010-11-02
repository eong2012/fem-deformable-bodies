#include "windows.h"
#define GLUT_DISABLE_ATEXIT_HACK
#define GLEW_STATIC

#include <GL/glew.h>
#include <GL/glut.h>

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

    void AddTetrahedron(vector<Vector3<float> > vertices);
    void AddHalfFace(unsigned int faceIndex,unsigned int vertexIndex);
    void AddHalfEdgePair(unsigned int vertexIndex1,unsigned int vertexIndex2,unsigned int &edgeIndex1,unsigned int &edgeIndex2);
    bool AddVertex(Vector3<float> vertexPos, unsigned int &index);
    void AddFace(unsigned int vertexIndex1,unsigned int vertexIndex2, unsigned int vertexIndex3,unsigned int &faceIndex);
    Vector3<float> FaceNormal(unsigned int faceIndex);
    void setTetraGeometry(vector<Vector3<float> > vertices, unsigned int &index1,unsigned int &index2,unsigned int &index3,unsigned int &index4);
    void Render();

private:

    vector<Vertex> mVertices;
    vector<HalfEdge> mHalfEdges;
    vector<Face> mFaces;
    vector<HalfFace> mHalfFaces;



};
