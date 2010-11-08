#include "TetrahedMesh.h"
#include <ctime>
#include <cstdlib>
using namespace std;
TetrahedMesh::TetrahedMesh()
{}
TetrahedMesh::~TetrahedMesh(){}

void TetrahedMesh::AddTetrahedron(vector<Vector3<float> > vertices)
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

	//Add a tetrahed which contains pointers to all spanning faces
	Tetrahed temp;
	temp.setFaceInd1(faceIndex1);
	temp.setFaceInd2(faceIndex2);
	temp.setFaceInd3(faceIndex3);
	temp.setFaceInd4(faceIndex4);
	mTetraheds.push_back(temp);


}

void TetrahedMesh::setTetraGeometry(vector<Vector3<float> > vertices, unsigned int &index1,unsigned int &index2,unsigned int &index3,unsigned int &index4)
{
//    unsigned int index1,index2,index3,index4;
    bool success = true;
    success &= AddVertex(vertices[0],index1);
    success &= AddVertex(vertices[1],index2);
    success &= AddVertex(vertices[2],index3);
    success &= AddVertex(vertices[3],index4);


}
void TetrahedMesh::AddFace(unsigned int vertexIndex1,unsigned int vertexIndex2, unsigned int vertexIndex3,unsigned int &faceIndex)
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

    vector<Face>::iterator iter;
	iter = mFaces.begin();
	int count = 0;
	while(iter != mFaces.end()){

    HalfEdge* edge = &mHalfEdges[iter->getEdgeInd()];

    Vertex & v1 = mVertices[edge->getVertexInd()];
    edge = &mHalfEdges[edge->getNextInd()];

    Vertex & v2 = mVertices[edge->getVertexInd()];
    edge = &mHalfEdges[edge->getNextInd()];

    Vertex & v3 = mVertices[edge->getVertexInd()];


	if (mVertices[vertexIndex1].getPosition() == v1.getPosition() || mVertices[vertexIndex1].getPosition() == v2.getPosition() || mVertices[vertexIndex1].getPosition() == v3.getPosition()){
		if(mVertices[vertexIndex2].getPosition() == v1.getPosition() || mVertices[vertexIndex2].getPosition() == v2.getPosition() || mVertices[vertexIndex2].getPosition() == v3.getPosition()){
			if(mVertices[vertexIndex3].getPosition() == v1.getPosition() || mVertices[vertexIndex3].getPosition() == v2.getPosition() || mVertices[vertexIndex3].getPosition() == v3.getPosition()){

				face.setOppositeFaceInd(count);
				iter->setOppositeFaceInd(mFaces.size());
				//cout <<  "Count: " << count << endl;
				break;
			}

		}
		}
		count++;

		iter += 1;

	}

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


void TetrahedMesh::AddHalfEdgePair(unsigned int vertexIndex1,unsigned int vertexIndex2,unsigned int &edgeIndex1,unsigned int &edgeIndex2)
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
bool TetrahedMesh::AddVertex(Vector3<float> vertexPos, unsigned int &index)
{

	vector<Vertex>::iterator iter;
	iter = mVertices.begin();
	int count = 0;
	while(iter != mVertices.end()){

		if (iter->getPosition() == vertexPos) {

			index = count;

			return true;
		}
		count++;
		iter += 1;
	}

  Vertex vert;
  vert.setPosition(vertexPos);
  mVertices.push_back(vert); // add it to the vertex list
  index = mVertices.size()-1;
  return true;
}

Vector3<float> TetrahedMesh::FaceNormal(unsigned int faceIndex)
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

void TetrahedMesh::Render(int mode)
{

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

		if(mode == 2){
			if (face.getOppositeFaceInd() == -1) {
			glColor3f(1.0f, 1.0f, 1.0f);

			glNormal3f(face.getNormal()[0], face.getNormal()[1], face.getNormal()[2]);

			glBegin(GL_TRIANGLES);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glEnd();
		}
		}

		if(mode == 1){
			if (face.getOppositeFaceInd() != -1) {
			glColor3f(1.0f, 1.0f, 1.0f);

			glNormal3f(face.getNormal()[0], face.getNormal()[1], face.getNormal()[2]);

			glBegin(GL_TRIANGLES);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glEnd();
			}
		}

		if(mode == 0){

			glColor3f(1.0f, 1.0f, 1.0f);

			glNormal3f(face.getNormal()[0], face.getNormal()[1], face.getNormal()[2]);

			glBegin(GL_TRIANGLES);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glEnd();

		}

  }



}

void TetrahedMesh::RenderNormals(int mode) {


	const int numTriangles = mFaces.size();

    for (int i = 0; i < numTriangles; i++){

    Face & face = mFaces[i];

    HalfEdge* edge = &mHalfEdges[face.getEdgeInd()];

    Vertex & v1 = mVertices[edge->getVertexInd()];
    edge = &mHalfEdges[edge->getNextInd()];

    Vertex & v2 = mVertices[edge->getVertexInd()];
    edge = &mHalfEdges[edge->getNextInd()];

    Vertex & v3 = mVertices[edge->getVertexInd()];

	Vector3<float> vCenter, normal;
	normal = face.getNormal();
	normal.Normalize();
	vCenter = (v1.getPosition()+v2.getPosition()+v3.getPosition())/3.0;

		if (mode == 0) {
			if (face.getOppositeFaceInd() != -1) {
			glColor3f(1.0f, 0.0f, 1.0f);
			glBegin(GL_LINES);
			glVertex3f(vCenter[0], vCenter[1], vCenter[2]);
			glVertex3f(vCenter[0]+normal[0]/4, vCenter[1]+normal[1]/4, vCenter[2]+normal[2]/4);

			glEnd();
			}

		}

		if (mode == 1) {
			if (face.getOppositeFaceInd() == -1) {
			glColor3f(1.0f, 0.0f, 1.0f);
			glBegin(GL_LINES);
			glVertex3f(vCenter[0], vCenter[1], vCenter[2]);
			glVertex3f(vCenter[0]+normal[0]/4, vCenter[1]+normal[1]/4, vCenter[2]+normal[2]/4);

			glEnd();
			}

		}

		if (mode == 2) {

			glColor3f(1.0f, 0.0f, 1.0f);
			glBegin(GL_LINES);
			glVertex3f(vCenter[0], vCenter[1], vCenter[2]);
			glVertex3f(vCenter[0]+normal[0]/4, vCenter[1]+normal[1]/4, vCenter[2]+normal[2]/4);

			glEnd();


		}
  }

}

void TetrahedMesh::RenderEdges(int mode) {

	const int numTriangles = mFaces.size();

    for (int i = 0; i < numTriangles; i++){

    Face & face = mFaces[i];


    HalfEdge* edge = &mHalfEdges[face.getEdgeInd()];

    Vertex & v1 = mVertices[edge->getVertexInd()];
    edge = &mHalfEdges[edge->getNextInd()];

    Vertex & v2 = mVertices[edge->getVertexInd()];
    edge = &mHalfEdges[edge->getNextInd()];

    Vertex & v3 = mVertices[edge->getVertexInd()];

		if (mode == 0) {
			if (face.getOppositeFaceInd() != -1) {
			glColor3f(0.2f, 0.0f, 1.0f);

			glBegin(GL_LINES);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glEnd();
			}

		}

		if (mode == 1) {
			if (face.getOppositeFaceInd() == -1) {
			glColor3f(0.2f, 0.0f, 1.0f);

			glBegin(GL_LINES);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glEnd();
			}

		}

		if (mode == 2) {

			glColor3f(0.2f, 0.0f, 1.0f);

			glBegin(GL_LINES);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v2.getPosition()[0], v2.getPosition()[1], v2.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glVertex3f(v3.getPosition()[0], v3.getPosition()[1], v3.getPosition()[2]);
			glVertex3f(v1.getPosition()[0], v1.getPosition()[1], v1.getPosition()[2]);
			glEnd();


		}
	}



}

//Get the vertex array for use to create a texture
float* TetrahedMesh::GetVertexArray()
{
    float *vertexArray = new float[4*mVertices.size()];

    for(int i = 0; i < mVertices.size(); i++)
    {
        vertexArray[i*4] = mVertices[i].getPosition()[0];
        vertexArray[i*4+1] = mVertices[i].getPosition()[1];
        vertexArray[i*4+2] = mVertices[i].getPosition()[2];
        vertexArray[i*4+3] = 1.0f;
    }
    return vertexArray;
}

int TetrahedMesh::GetVertexArraySize()
{
    return 4*mVertices.size();
}


vector<arma::Mat<double> > TetrahedMesh::getVertexPosition(unsigned int tetraIndex){
    Tetrahed tempTetra;
    tempTetra = mTetraheds.at(tetraIndex);

    vector<unsigned int> faceIndices = tempTetra.getFaceInd();
    set<unsigned int> vertexIndices;

	vector<arma::Mat<double> > ret;
    //vector< Vector3<float> > ret;

    for(int i=0; i<4;i++){

        Face & face = mFaces[faceIndices[i]];

        HalfEdge* edge = &mHalfEdges[face.getEdgeInd()];

        vertexIndices.insert(edge->getVertexInd());

        edge = &mHalfEdges[edge->getNextInd()];
        vertexIndices.insert(edge->getVertexInd());

        edge = &mHalfEdges[edge->getNextInd()];
        vertexIndices.insert(edge->getVertexInd());
    }

    set<unsigned int>::iterator iter = vertexIndices.begin();

    while(iter!=vertexIndices.end())
    {
		arma::Mat<double> temp(4,1);
		temp(0) = mVertices[*iter].getPosition()[0];
		temp(1) = mVertices[*iter].getPosition()[1];
		temp(2) = mVertices[*iter].getPosition()[2];
		temp(3) = (*iter);
        ret.push_back(temp);
		//cout << *iter << endl;
        iter++;
    }

    return ret;
}

int TetrahedMesh::getNrOfTetrahedra(){
    return mTetraheds.size();
}

void TetrahedMesh::updatePosition(unsigned int tetraIndex, vector<Vector3<float> > vertexPos)
{
    Tetrahed tempTetra;
    tempTetra = mTetraheds.at(tetraIndex);

    vector<unsigned int> faceIndices = tempTetra.getFaceInd();
    set<unsigned int> vertexIndices;

     for(int i=0; i<4;i++){

        Face & face = mFaces[faceIndices[i]];

        HalfEdge* edge = &mHalfEdges[face.getEdgeInd()];

        vertexIndices.insert(edge->getVertexInd());

        edge = &mHalfEdges[edge->getNextInd()];
        vertexIndices.insert(edge->getVertexInd());

        edge = &mHalfEdges[edge->getNextInd()];
        vertexIndices.insert(edge->getVertexInd());
    }

    set<unsigned int>::iterator iter = vertexIndices.begin();
    int count = 0;
    while(iter!=vertexIndices.end())
    {
        mVertices[*iter].setPosition(vertexPos[count]);
        iter++;
        count++;
    }

}

void TetrahedMesh::printVertices(int ind)
{
    cout << mVertices[ind].getPosition() << endl;
}
