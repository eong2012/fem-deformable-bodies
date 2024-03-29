#include "HalfFaceMesh.h"
#include <ctime>
#include <cstdlib>
using namespace std;
HalfFaceMesh::HalfFaceMesh()
{}
HalfFaceMesh::~HalfFaceMesh(){}

void HalfFaceMesh::AddTetrahedron(vector<Vector3<float> > vertices)
{

    //set the vertex index
    unsigned int vertexIndex1, vertexIndex2, vertexIndex3, vertexIndex4;
    setTetraGeometry(vertices,vertexIndex1, vertexIndex2, vertexIndex3, vertexIndex4);

    //Add each tetrahedron face, send in vertex indices counter-clockwise
    unsigned int faceIndex1, faceIndex2,faceIndex3,faceIndex4;

	AddFace(vertexIndex1, vertexIndex4, vertexIndex3,faceIndex1);
    AddFace(vertexIndex1, vertexIndex3, vertexIndex2,faceIndex2);
    AddFace(vertexIndex2, vertexIndex3, vertexIndex4,faceIndex3);
    AddFace(vertexIndex1, vertexIndex2, vertexIndex4,faceIndex4);

    //AddHalfFace
    //AddHalfFace(faceIndex1, vertexIndex4);

}

void HalfFaceMesh::setTetraGeometry(vector<Vector3<float> > vertices, unsigned int &index1,unsigned int &index2,unsigned int &index3,unsigned int &index4)
{
//    unsigned int index1,index2,index3,index4;
    bool success = true;
    success &= AddVertex(vertices[0],index1);
    success &= AddVertex(vertices[1],index2);
    success &= AddVertex(vertices[2],index3);
    success &= AddVertex(vertices[3],index4);


}
void HalfFaceMesh::AddFace(unsigned int vertexIndex1,unsigned int vertexIndex2, unsigned int vertexIndex3,unsigned int &faceIndex)
{
    //add halfedges
    unsigned int edgeIndex1,edgeIndex2,edgeIndex3,edgeIndex4,edgeIndex5,edgeIndex6;
    AddHalfEdgePair(vertexIndex1,vertexIndex2,edgeIndex1,edgeIndex2);
    AddHalfEdgePair(vertexIndex2,vertexIndex3,edgeIndex3,edgeIndex4);
    AddHalfEdgePair(vertexIndex3,vertexIndex1,edgeIndex5,edgeIndex6);

    // Connect inner ring
    mHalfEdges[edgeIndex1].setNextInd(edgeIndex3);
    mHalfEdges[edgeIndex1].setPrevInd(edgeIndex5);

    mHalfEdges[edgeIndex3].setNextInd(edgeIndex5);
    mHalfEdges[edgeIndex3].setPrevInd(edgeIndex1);

    mHalfEdges[edgeIndex5].setNextInd(edgeIndex1);
    mHalfEdges[edgeIndex5].setPrevInd(edgeIndex3);

    // Finally, create the face, don't forget to set the normal

    Face face;

    // Connect the face to an edge
    face.setEdgeInd(edgeIndex1);

    mFaces.push_back(face);

    // Compute and assign a normal
    mFaces.back().setNormal(FaceNormal(mFaces.size() - 1));

    // All half-edges share the same left face (previously added)
    mHalfEdges[edgeIndex1].setFaceInd(mFaces.size() - 1);
    mHalfEdges[edgeIndex2].setFaceInd(mFaces.size() - 1);
    mHalfEdges[edgeIndex3].setFaceInd(mFaces.size() - 1);
    //also send back faceIndex
    faceIndex = mFaces.size() - 1;

}

void HalfFaceMesh::AddHalfFace(unsigned int faceIndex,unsigned int vertexIndex)
{
    HalfFace halfFace;
    //Set the halfface vertex
    halfFace.setVertexInd(vertexIndex);

    //set face
    halfFace.setFaceInd(faceIndex);

    //Set the right edge, the one going from the face to vertex4
    unsigned int i = mFaces[faceIndex].getEdgeInd();
    unsigned int j = mHalfEdges[i].getPairInd();
    unsigned int k = mHalfEdges[j].getNextInd();
    halfFace.setEdgeInd(k);

    //Set radial somehow :)


}

void HalfFaceMesh::AddHalfEdgePair(unsigned int vertexIndex1,unsigned int vertexIndex2,unsigned int &edgeIndex1,unsigned int &edgeIndex2)
{
  //Set index
  edgeIndex1 = mHalfEdges.size();
  edgeIndex2 = edgeIndex1+1;

  // Create edges and set pair index
  HalfEdge halfEdge1, halfEdge2;
  halfEdge1.setPairInd(edgeIndex2);
  halfEdge2.setPairInd(edgeIndex1);

  // Connect the edges to the verts
  halfEdge1.setVertexInd(vertexIndex1);
  halfEdge2.setVertexInd(vertexIndex2);

  // Connect the verts to the edges
  mVertices[vertexIndex1].setEdgeInd(edgeIndex1);
  mVertices[vertexIndex2].setEdgeInd(edgeIndex2);

   // Store the edges in mEdges
  mHalfEdges.push_back(halfEdge1);
  mHalfEdges.push_back(halfEdge2);

  // Store the first edge in the map as an OrderedPair
 // OrderedPair op(v1, v2);
}
bool HalfFaceMesh::AddVertex(Vector3<float> vertexPos, unsigned int &index)
{
  /*  std::map<Vector3<float>,unsigned int>::iterator it = mUniqueVerts.find(v);
  if (it != mUniqueVerts.end()){
    indx = (*it).second; // get the index of the already existing vertex
    return false;
  }

  mUniqueVerts[vertexPos] = indx = GetNumVerts(); // op. [ ] constructs a new entry in map*/

  Vertex vert;
  vert.setPosition(vertexPos);
  mVertices.push_back(vert); // add it to the vertex list
  index = mVertices.size()-1;
  return true;
}

Vector3<float> HalfFaceMesh::FaceNormal(unsigned int faceIndex)
{
  unsigned int edgeInd1 = mFaces[faceIndex].getEdgeInd();
  unsigned int vertInd1 = mHalfEdges[edgeInd1].getVertexInd();
  const Vector3<float> &p1 = mVertices[vertInd1].getPosition();

  unsigned int edgeInd2 = mHalfEdges[edgeInd1].getNextInd();
  unsigned int vertInd2 = mHalfEdges[edgeInd2].getVertexInd();
  const Vector3<float> &p2 = mVertices[vertInd2].getPosition();

  unsigned int edgeInd3 = mHalfEdges[edgeInd2].getNextInd();
  unsigned int vertInd3 = mHalfEdges[edgeInd3].getVertexInd();
  const Vector3<float> &p3 = mVertices[vertInd3].getPosition();


  const Vector3<float> e1 = p2-p1;
  const Vector3<float> e2 = p3-p1;
  return Cross(e1, e2).Normalize();
}

void HalfFaceMesh::Render()
{
     glDisable(GL_LIGHTING);
// Draw geometry

  const int numTriangles = mFaces.size();
  for (int i = 0; i < numTriangles; i++){

    Face & face = mFaces[i];


    HalfEdge* edge = &mHalfEdges[face.getEdgeInd()];

    Vertex & v1 = mVertices[edge->getVertexInd()];
    edge = &mHalfEdges[edge->getNextInd()];

    Vertex & v2 = mVertices[edge->getVertexInd()];
    edge = &mHalfEdges[edge->getNextInd()];

    Vertex & v3 = mVertices[edge->getVertexInd()];

    glColor3f(1.0f, 1.0f, 1.0f);
    glNormal3f(face.getNormal()[0], face.getNormal()[1], face.getNormal()[2]);

	glBegin(GL_TRIANGLES);
    glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
    glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
    glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
	glEnd();

	glColor3f(0.2f, 0.0f, 1.0f);
	/*
	glBegin(GL_LINES);
    glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
    glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
	glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
	glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
	glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
	glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
	glEnd();
	*/

	if(i == 0) {
	glColor3f(1.0f, 0.0f, 0.0f);
	} else if(i == 1) {
	glColor3f(0.0f, 1.0f, 0.0f);
	}
	 else if(i == 2) {
	glColor3f(0.0f, 0.0f, 1.0f);
	}
	  else if(i == 3) {
	glColor3f(1.0f, 1.0f, 0.0f);
	}

	Vector3<float> vCenter, normal;
	normal = face.getNormal();
	normal.Normalize();
	vCenter = (v1.getPosition()+v2.getPosition()+v3.getPosition())/3.0;

	glBegin(GL_LINES);
	glVertex3f(vCenter[0], vCenter[1], vCenter[2]);
	glVertex3f(vCenter[0]+normal[0]/4, vCenter[1]+normal[1]/4, vCenter[2]+normal[2]/4);

	glEnd();


  }



}
