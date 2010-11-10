#include "windows.h"
#define GLUT_DISABLE_ATEXIT_HACK
#define GLEW_STATIC

#include <gl/glew.h>
#include <gl/glut.h>

#include "Tetrahed.h"
#include "Face.h"
#include "Vertex.h"
#include "HalfEdge.h"
#include "armadillo"


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
	void buildTetrahedonMesh(vector<Vector3<float> > vertices, vector<Tetrahed> *mtetrahedsTemp, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp, bool subdivide);
    void AddHalfEdgePair(unsigned int vertexIndex1,unsigned int vertexIndex2,unsigned int &edgeIndex1,unsigned int &edgeIndex2, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp);
    bool AddVertex(Vector3<float> vertexPos, unsigned int &index, vector<Vertex> *mVerticesTemp);
    void AddFace(unsigned int vertexIndex1,unsigned int vertexIndex2, unsigned int vertexIndex3,unsigned int &faceIndex, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp);

    Vector3<float> FaceNormal(unsigned int faceIndex, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp);
    void setTetraGeometry(vector<Vector3<float> > vertices, unsigned int &index1,unsigned int &index2,unsigned int &index3,unsigned int &index4, vector<Vertex> *mVerticesTemp);
	

    void Render(int mode);
	void RenderNormals(int mode);
	void RenderEdges(int mode);

	vector< arma::Mat<double> > getVertexPosition(unsigned int tetraIndex);
	int getNrOfNodes() {return mVertices->size();}
	int getNrOfTetrahedra();

	float* GetVertexArray();
	int GetVertexArraySize();

	void subdivide();

	void updatePosition(unsigned int tetraIndex, vector<Vector3<float> > vertexPos);

	//temp debug thing
	void printVertices(int ind);
	vector<Vertex> *mVertices;

private:

    
    vector<HalfEdge> *mHalfEdges;
    vector<Face> *mFaces;
    vector<Tetrahed> *mTetraheds;

};

#endif
