//#include "windows.h"
#define GLUT_DISABLE_ATEXIT_HACK
#define GLEW_STATIC

#include <GL/glew.h>
#include <GL/glut.h>

#include "Tetrahed.h"
#include "Face.h"
#include "Vertex.h"
#include "HalfEdge.h"

#include <vector>
#include <set>

#ifndef TETRAHEDMESH_H
#define TETRAHEDMESH_H

using namespace std;

class TetrahedMesh{
public:
    TetrahedMesh();
    ~TetrahedMesh();

    void AddTetrahedron(vector<Vector3<float> > vertices);
    void AddHalfEdgePair(unsigned int vertexIndex1,unsigned int vertexIndex2,unsigned int &edgeIndex1,unsigned int &edgeIndex2);
    bool AddVertex(Vector3<float> vertexPos, unsigned int &index);
    void AddFace(unsigned int vertexIndex1,unsigned int vertexIndex2, unsigned int vertexIndex3,unsigned int &faceIndex);
    Vector3<float> FaceNormal(unsigned int faceIndex);
    void setTetraGeometry(vector<Vector3<float> > vertices, unsigned int &index1,unsigned int &index2,unsigned int &index3,unsigned int &index4);
    void Render(int mode);
	void RenderNormals(int mode);
	void RenderEdges(int mode);

	vector< Vector3<float> > getVertexPosition(unsigned int tetraIndex);
	int getNrOfTetrahedra();

	float* GetVertexArray();
	int GetVertexArraySize();

	void updatePosition(unsigned int tetraIndex, vector<Vector3<float> > vertexPos);

	//temp debug thing
	void printVertices(int ind);

private:

    vector<Vertex> mVertices;
    vector<HalfEdge> mHalfEdges;
    vector<Face> mFaces;
    vector<Tetrahed> mTetraheds;

};

#endif