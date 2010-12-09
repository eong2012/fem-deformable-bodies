
#define GLUT_DISABLE_ATEXIT_HACK
#define GLEW_STATIC

#include <GL/glew.h>
#include <GL/glut.h>

#include "Tetrahed.h"
#include "Face.h"
#include "Vertex.h"
#include "HalfEdge.h"
#include "../../usr/include/armadillo.h"


#include <vector>
#include <set>

#ifndef TETRAHEDMESH_H
#define TETRAHEDMESH_H

using namespace std;

class TetrahedMesh{
public:
    TetrahedMesh();
    ~TetrahedMesh();

    void AddTetrahedron(vector<arma::Mat<double> > vertices);
	void buildTetrahedonMesh(vector<arma::Mat<double> > vertices, vector<Tetrahed> *mtetrahedsTemp, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp, bool subdivide);
    void AddHalfEdgePair(unsigned int vertexIndex1,unsigned int vertexIndex2,unsigned int &edgeIndex1,unsigned int &edgeIndex2, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp);
    bool AddVertex(arma::Mat<double> vertexPos, unsigned int &index, vector<Vertex> *mVerticesTemp);
    void AddFace(unsigned int vertexIndex1,unsigned int vertexIndex2, unsigned int vertexIndex3,unsigned int &faceIndex, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp);


    arma::Mat<double> FaceNormal(unsigned int faceIndex, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp);
    void setTetraGeometry(vector<arma::Mat<double> > vertices, unsigned int &index1,unsigned int &index2,unsigned int &index3,unsigned int &index4, vector<Vertex> *mVerticesTemp);

	void crackStructure(unsigned int ind, arma::Mat<double> eigVec);
	vector<unsigned int> getAdjecentTetraheds(unsigned int ind, arma::Mat<double> v);
	void getAdjecent(unsigned int currentind, unsigned int baseind, vector<unsigned int> *tetraIndices, arma::Mat<double> v, int count);
	arma::Mat<double> tetCenterOfMass(unsigned int tInd);
	arma::Mat<double> determineLocationOfCrack(unsigned int tInd, vector<unsigned int> adjecentList, arma::Mat<double> v, arma::Mat<double> eigVec);
	void deconnect(vector<unsigned int> adjecentList, arma::Mat<double> v);

    void Render(int mode);
	void RenderNormals(int mode);
	void RenderEdges(int mode);

	vector< arma::Mat<double> > getVertexPosition(unsigned int tetraIndex);
	int getNrOfNodes() {return mVertices->size();}
	int getNrOfTetrahedra();

	bool vecEquals(arma::Mat<double> A, arma::Mat<double> B);

	float* GetVertexArray();
	int GetVertexArraySize();

	void subdivide();

	void updatePosition(unsigned int tetraIndex, vector<arma::Mat<double> > vertexPos);

	arma::Mat<double> pickNode();
	void pickNextNode();
	unsigned int getCurrentNode();

	//temp debug thing
	void printVertices(int ind);
	vector<Vertex> *mVertices;
	vector<unsigned int> *mVertexIndexOrder;
	vector<Vertex> mMaterialCoords;

private:


    vector<HalfEdge> *mHalfEdges;
    vector<Face> *mFaces;
    vector<Tetrahed> *mTetraheds;

	unsigned int currentNode;

};

#endif
