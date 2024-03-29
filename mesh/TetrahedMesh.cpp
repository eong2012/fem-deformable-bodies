#include "TetrahedMesh.h"
#include <ctime>
#include <cstdlib>
using namespace std;
TetrahedMesh::TetrahedMesh()
{
	mTetraheds = new vector<Tetrahed>();
	mHalfEdges = new vector<HalfEdge>();
	mFaces = new vector<Face>();
	mVertices = new vector<Vertex>();
	mVertexIndexOrder = new vector<unsigned int>();
	this->currentNode = 0;
	srand(time(0));



}
TetrahedMesh::~TetrahedMesh(){

    delete mTetraheds;
    delete mHalfEdges;
    delete mFaces;
    delete mVertices;
    delete mVertexIndexOrder;
}

void TetrahedMesh::AddTetrahedron(vector<arma::Mat<double> > vertices ) {

	buildTetrahedonMesh(vertices, mTetraheds, mFaces, mHalfEdges, mVertices, false);
}

void TetrahedMesh::buildTetrahedonMesh(vector<arma::Mat<double> > vertices, vector<Tetrahed> *mtetrahedsTemp, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp, bool subdivide)
{
    //set the vertex index
    unsigned int vertexIndex1, vertexIndex2, vertexIndex3, vertexIndex4;
    setTetraGeometry(vertices,vertexIndex1, vertexIndex2, vertexIndex3, vertexIndex4, mVerticesTemp);

    mVertexIndexOrder->push_back(vertexIndex1);
    mVertexIndexOrder->push_back(vertexIndex2);
    mVertexIndexOrder->push_back(vertexIndex3);
    mVertexIndexOrder->push_back(vertexIndex4);

    //Add each tetrahedron face, send in vertex indices counter-clockwise
    unsigned int faceIndex1, faceIndex2,faceIndex3,faceIndex4;

	AddFace(vertexIndex1, vertexIndex4, vertexIndex3,faceIndex1, mFacesTemp, mEdgesTemp, mVerticesTemp);
    AddFace(vertexIndex1, vertexIndex3, vertexIndex2,faceIndex2, mFacesTemp, mEdgesTemp, mVerticesTemp);
    AddFace(vertexIndex2, vertexIndex3, vertexIndex4,faceIndex3, mFacesTemp, mEdgesTemp, mVerticesTemp);
    AddFace(vertexIndex1, vertexIndex2, vertexIndex4,faceIndex4, mFacesTemp, mEdgesTemp, mVerticesTemp);



	//Add a tetrahed which contains pointers to all spanning faces
	Tetrahed* temp = new Tetrahed();
	temp->setFaceInd1(faceIndex1);
	temp->setFaceInd2(faceIndex2);
	temp->setFaceInd3(faceIndex3);
	temp->setFaceInd4(faceIndex4);






	mtetrahedsTemp->push_back((*temp));
	int sis = mtetrahedsTemp->size();
	mFacesTemp->at(faceIndex1).setTetrahedInd(mtetrahedsTemp->size()-1);
	mFacesTemp->at(faceIndex2).setTetrahedInd(mtetrahedsTemp->size()-1);
	mFacesTemp->at(faceIndex3).setTetrahedInd(mtetrahedsTemp->size()-1);
	mFacesTemp->at(faceIndex4).setTetrahedInd(mtetrahedsTemp->size()-1);
	delete temp;

	if (subdivide == false){
	mTetraheds = mtetrahedsTemp;
	mHalfEdges = mEdgesTemp;
	mVertices = mVerticesTemp;
	mFaces = mFacesTemp;
	}

}

void TetrahedMesh::setTetraGeometry(vector<arma::Mat<double> > vertices, unsigned int &index1,unsigned int &index2,unsigned int &index3,unsigned int &index4, vector<Vertex> *mVertices)
{
//    unsigned int index1,index2,index3,index4;
    bool success = true;
    success &= AddVertex(vertices[0],index1, mVertices);
    success &= AddVertex(vertices[1],index2, mVertices);
    success &= AddVertex(vertices[2],index3, mVertices);
    success &= AddVertex(vertices[3],index4, mVertices);


}
void TetrahedMesh::AddFace(unsigned int vertexIndex1,unsigned int vertexIndex2, unsigned int vertexIndex3,unsigned int &faceIndex, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp)
{
    //add halfedges
    unsigned int edgeIndex1,edgeIndex2,edgeIndex3,edgeIndex4,edgeIndex5,edgeIndex6;
    AddHalfEdgePair(vertexIndex1,vertexIndex2,edgeIndex1,edgeIndex2, mEdgesTemp,mVerticesTemp);
    AddHalfEdgePair(vertexIndex2,vertexIndex3,edgeIndex3,edgeIndex4, mEdgesTemp,mVerticesTemp);
    AddHalfEdgePair(vertexIndex3,vertexIndex1,edgeIndex5,edgeIndex6, mEdgesTemp,mVerticesTemp);

    // Connect inner ring
    mEdgesTemp->at(edgeIndex1).setNextInd(edgeIndex3);
    mEdgesTemp->at(edgeIndex1).setPrevInd(edgeIndex5);

    mEdgesTemp->at(edgeIndex3).setNextInd(edgeIndex5);
    mEdgesTemp->at(edgeIndex3).setPrevInd(edgeIndex1);

    mEdgesTemp->at(edgeIndex5).setNextInd(edgeIndex1);
    mEdgesTemp->at(edgeIndex5).setPrevInd(edgeIndex3);

  	 // Finally, create the face, don't forget to set the normal
	Face *face = new Face();


    vector<Face>::iterator iter;
	iter = mFacesTemp->begin();
	int count = 0;
	while(iter != mFacesTemp->end()){

		HalfEdge* edge = &mEdgesTemp->at(iter->getEdgeInd());

		Vertex & v1 = mVerticesTemp->at(edge->getVertexInd());
		edge = &mEdgesTemp->at(edge->getNextInd());

		Vertex & v2 = mVerticesTemp->at(edge->getVertexInd());
		edge = &mEdgesTemp->at(edge->getNextInd());

		Vertex & v3 = mVerticesTemp->at(edge->getVertexInd());





	if (vecEquals(mVerticesTemp->at(vertexIndex1).getPosition(), v1.getPosition()) || vecEquals(mVerticesTemp->at(vertexIndex1).getPosition(), v2.getPosition()) || vecEquals(mVerticesTemp->at(vertexIndex1).getPosition(), v3.getPosition())){
		if(vecEquals(mVerticesTemp->at(vertexIndex2).getPosition(),v1.getPosition()) || vecEquals(mVerticesTemp->at(vertexIndex2).getPosition(), v2.getPosition()) || vecEquals(mVerticesTemp->at(vertexIndex2).getPosition(),v3.getPosition())){
			if(vecEquals(mVerticesTemp->at(vertexIndex3).getPosition(), v1.getPosition()) || vecEquals(mVerticesTemp->at(vertexIndex3).getPosition(),v2.getPosition()) || vecEquals(mVerticesTemp->at(vertexIndex3).getPosition(), v3.getPosition())){

				face->setOppositeFaceInd(count);
				iter->setOppositeFaceInd(mFacesTemp->size());
				//cout <<  "Count: " << count << endl;
				break;
			}

		}
		}
		count++;

		iter += 1;

	}

    // Connect the face to an edge
    face->setEdgeInd(edgeIndex1);

    mFacesTemp->push_back((*face));
	delete face;

    // Compute and assign a normal
    mFacesTemp->back().setNormal(FaceNormal(mFacesTemp->size() - 1, mFacesTemp, mEdgesTemp, mVerticesTemp));

    // All half-edges share the same left face (previously added)
    mEdgesTemp->at(edgeIndex1).setFaceInd(mFacesTemp->size() - 1);
    mEdgesTemp->at(edgeIndex2).setFaceInd(mFacesTemp->size() - 1);
    mEdgesTemp->at(edgeIndex3).setFaceInd(mFacesTemp->size() - 1);
    //also send back faceIndex
    faceIndex = mFacesTemp->size() - 1;



}


void TetrahedMesh::AddHalfEdgePair(unsigned int vertexIndex1,unsigned int vertexIndex2,unsigned int &edgeIndex1,unsigned int &edgeIndex2, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp)
{
  //Set index
  edgeIndex1 = mEdgesTemp->size();
  edgeIndex2 = edgeIndex1+1;

  // Create edges and set pair index
  HalfEdge *halfEdge1, *halfEdge2;

  halfEdge1 = new HalfEdge();
  halfEdge2 = new HalfEdge();

  halfEdge1->setPairInd(edgeIndex2);
  halfEdge2->setPairInd(edgeIndex1);

  // Connect the edges to the verts
  halfEdge1->setVertexInd(vertexIndex1);
  halfEdge2->setVertexInd(vertexIndex2);

  // Connect the verts to the edges
  mVerticesTemp->at(vertexIndex1).setEdgeInd(edgeIndex1);
  mVerticesTemp->at(vertexIndex2).setEdgeInd(edgeIndex2);

   // Store the edges in mEdges
  mEdgesTemp->push_back(*halfEdge1);
  mEdgesTemp->push_back(*halfEdge2);
  delete halfEdge1;
  delete halfEdge2;

  // Store the first edge in the map as an OrderedPair
 // OrderedPair op(v1, v2);
}
bool TetrahedMesh::AddVertex(arma::Mat<double> vertexPos, unsigned int &index, vector<Vertex> *mVerticesTemp)
{

	vector<Vertex>::iterator iter;
	iter = mVerticesTemp->begin();
	int count = 0;
	while(iter != mVerticesTemp->end()){

		if (vecEquals(iter->getPosition(), vertexPos)) {

			index = count;

			return true;
		}
		count++;
		iter += 1;
	}

  Vertex vert;
  vert.setPosition(vertexPos);
  mVerticesTemp->push_back(vert); // add it to the vertex list
  index = mVerticesTemp->size()-1;
  return true;
}

arma::Mat<double> TetrahedMesh::FaceNormal(unsigned int faceIndex, vector<Face> *mFacesTemp, vector<HalfEdge> *mEdgesTemp, vector<Vertex> *mVerticesTemp)
{
  unsigned int edgeInd1 = mFacesTemp->at(faceIndex).getEdgeInd();
  unsigned int vertInd1 = mEdgesTemp->at(edgeInd1).getVertexInd();
  const arma::Mat<double> &p1 = mVerticesTemp->at(vertInd1).getPosition();

  unsigned int edgeInd2 = mEdgesTemp->at(edgeInd1).getNextInd();
  unsigned int vertInd2 = mEdgesTemp->at(edgeInd2).getVertexInd();
  const arma::Mat<double> &p2 = mVerticesTemp->at(vertInd2).getPosition();

  unsigned int edgeInd3 = mEdgesTemp->at(edgeInd2).getNextInd();
  unsigned int vertInd3 = mEdgesTemp->at(edgeInd3).getVertexInd();
  const arma::Mat<double> &p3 = mVerticesTemp->at(vertInd3).getPosition();


  const arma::Mat<double> e1 = p2-p1;
  const arma::Mat<double> e2 = p3-p1;
  arma::Mat<double> tmp;
  tmp = cross(e1, e2);
  tmp = tmp / norm(tmp,1);
  return tmp;
}

void TetrahedMesh::Render(int mode)
{

// Draw geometry

  const int numTriangles = mFaces->size();
  for (int i = 0; i < numTriangles; i++){

    Face & face = mFaces->at(i);
    face.setNormal(FaceNormal(i,mFaces, mHalfEdges, mVertices));
    HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

    Vertex & v1 = mVertices->at(edge->getVertexInd());
    edge = &mHalfEdges->at(edge->getNextInd());

    Vertex & v2 = mVertices->at(edge->getVertexInd());
    edge = &mHalfEdges->at(edge->getNextInd());

    Vertex & v3 = mVertices->at(edge->getVertexInd());

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


	const int numTriangles = mFaces->size();

    for (int i = 0; i < numTriangles; i++){

    Face & face = mFaces->at(i);

    HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

    Vertex & v1 = mVertices->at(edge->getVertexInd());
    edge = &mHalfEdges->at(edge->getNextInd());

    Vertex & v2 = mVertices->at(edge->getVertexInd());
    edge = &mHalfEdges->at(edge->getNextInd());

    Vertex & v3 = mVertices->at(edge->getVertexInd());

	arma::Mat<double> vCenter, normal;
	normal = face.getNormal();
	float scale= 80;
	vCenter = (v1.getPosition()+v2.getPosition()+v3.getPosition())/3.0;

		if (mode == 0) {
			if (face.getOppositeFaceInd() != -1) {
			glColor3f(1.0f, 0.0f, 1.0f);
			glBegin(GL_LINES);
			glVertex3f(vCenter[0], vCenter[1], vCenter[2]);
			glVertex3f(vCenter[0]+normal[0]/scale, vCenter[1]+normal[1]/scale, vCenter[2]+normal[2]/scale);

			glEnd();
			}

		}

		if (mode == 1) {
			if (face.getOppositeFaceInd() == -1) {
			glColor3f(1.0f, 0.0f, 1.0f);
			glBegin(GL_LINES);
			glVertex3f(vCenter[0], vCenter[1], vCenter[2]);
			glVertex3f(vCenter[0]+normal[0]/scale, vCenter[1]+normal[1]/scale, vCenter[2]+normal[2]/scale);

			glEnd();
			}

		}

		if (mode == 2) {

			glColor3f(1.0f, 0.0f, 1.0f);
			glBegin(GL_LINES);
			glVertex3f(vCenter[0], vCenter[1], vCenter[2]);
			glVertex3f(vCenter[0]+normal[0]/scale, vCenter[1]+normal[1]/scale, vCenter[2]+normal[2]/scale);

			glEnd();


		}
  }

}

void TetrahedMesh::RenderEdges(int mode) {

	const int numTriangles = mFaces->size();

    for (int i = 0; i < numTriangles; i++){

		Face & face = mFaces->at(i);


    HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

    Vertex & v1 = mVertices->at(edge->getVertexInd());
    edge = &mHalfEdges->at(edge->getNextInd());

    Vertex & v2 = mVertices->at(edge->getVertexInd());
    edge = &mHalfEdges->at(edge->getNextInd());

    Vertex & v3 = mVertices->at(edge->getVertexInd());

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
    float *vertexArray = new float[4*mVertices->size()];

    for(unsigned int i = 0; i < mVertices->size(); i++)
    {
        vertexArray[i*4] = mVertices->at(i).getPosition()[0];
        vertexArray[i*4+1] = mVertices->at(i).getPosition()[1];
        vertexArray[i*4+2] = mVertices->at(i).getPosition()[2];
        vertexArray[i*4+3] = 1.0f;
    }
    return vertexArray;
}

int TetrahedMesh::GetVertexArraySize()
{
    return 4*mVertices->size();
}


vector<arma::Mat<double> > TetrahedMesh::getVertexPosition(unsigned int tetraIndex){
    vector<arma::Mat<double> > ret;
    unsigned int vertexIndex[4];


    for(int i=0; i<4;i++){
        arma::Mat<double> tempVertexPos(4,1);
        vertexIndex[i] = mVertexIndexOrder->at(tetraIndex*4+i);

        tempVertexPos(0) = mVertices->at(vertexIndex[i]).getPosition()[0];
        tempVertexPos(1) = mVertices->at(vertexIndex[i]).getPosition()[1];
        tempVertexPos(2) = mVertices->at(vertexIndex[i]).getPosition()[2];
        tempVertexPos(3) = vertexIndex[i];
        ret.push_back(tempVertexPos);
    }
    return ret;


   /* Tetrahed tempTetra;
    tempTetra = mTetraheds->at(tetraIndex);

    vector<unsigned int> faceIndices = tempTetra.getFaceInd();
    set<unsigned int> vertexIndices;

	vector<arma::Mat<double> > ret;
    //vector< arma::Mat<double> > ret;

    for(int i=0; i<4;i++){

		Face & face = mFaces->at(faceIndices[i]);

        HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

        vertexIndices.insert(edge->getVertexInd());

        edge = &mHalfEdges->at(edge->getNextInd());
        vertexIndices.insert(edge->getVertexInd());

        edge = &mHalfEdges->at(edge->getNextInd());
        vertexIndices.insert(edge->getVertexInd());
    }

    set<unsigned int>::iterator iter = vertexIndices.begin();

    while(iter!=vertexIndices.end())
    {
		arma::Mat<double> temp(4,1);
		temp(0) = mVertices->at(*iter).getPosition()[0];
		temp(1) = mVertices->at(*iter).getPosition()[1];
		temp(2) = mVertices->at(*iter).getPosition()[2];
		temp(3) = (*iter);
        ret.push_back(temp);
		//cout << *iter << endl;
        iter++;
    }
*/

}

int TetrahedMesh::getNrOfTetrahedra(){
    return mTetraheds->size();
}

void TetrahedMesh::updatePosition(unsigned int tetraIndex, vector<arma::Mat<double> > vertexPos)
{
    Tetrahed tempTetra;
    tempTetra = mTetraheds->at(tetraIndex);

    vector<unsigned int> faceIndices = tempTetra.getFaceInd();
    set<unsigned int> vertexIndices;

     for(int i=0; i<4;i++){

        Face & face = mFaces->at(faceIndices[i]);

        HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

        vertexIndices.insert(edge->getVertexInd());

        edge = &mHalfEdges->at(edge->getNextInd());
        vertexIndices.insert(edge->getVertexInd());

        edge = &mHalfEdges->at(edge->getNextInd());
        vertexIndices.insert(edge->getVertexInd());
    }

    set<unsigned int>::iterator iter = vertexIndices.begin();
    int count = 0;
    while(iter!=vertexIndices.end())
    {
        mVertices->at(*iter).setPosition(vertexPos[count]);
        iter++;
        count++;
    }

}

void TetrahedMesh::printVertices(int ind)
{
    cout << mVertices->at(ind).getPosition() << endl;
}


bool TetrahedMesh::vecEquals(arma::Mat<double> A, arma::Mat<double> B){

      arma::umat C = (A == B);
      for(unsigned int i= 0; i<C.n_elem;i++){
          if(C(i)== 0) return false;
      }
      return true;

}


void TetrahedMesh::subdivide() {


	vector<HalfEdge> *newMHalfEdges = new vector<HalfEdge>();
    vector<Face> *newMFaces = new vector<Face>();
    vector<Tetrahed> *newMTetraheds = new vector<Tetrahed>();
	vector<Vertex> *newMVertices = new vector<Vertex>();

	for(unsigned int i = 0; i < mTetraheds->size(); i++) {

		arma::Mat<double> v5;
		v5 = arma::zeros(3,1);

		vector<arma::Mat<double> > verticesList;
		vector<arma::Mat<double> > verticesT = getVertexPosition(i);
		for (unsigned int j= 0; j < verticesT.size(); j++) {


			v5 += verticesT[j].rows(0,2);
			arma::Mat<double> temp;
			temp << verticesT[j](0) << verticesT[j](1) << verticesT[j](2);
			verticesList.push_back(temp);
		}

		v5 = v5/verticesT.size();
		verticesT.clear();
		arma::Mat<double> v5p;
		v5p << v5(0) << v5(1) << v5(2);
		vector<arma::Mat<double> > verticesList2;

		//TETRAHED 1

		if(i == 1 || i == 2) {
		verticesList2.push_back(verticesList.at(2));
		verticesList2.push_back(verticesList.at(1));
		verticesList2.push_back(verticesList.at(0));
		verticesList2.push_back(v5p);

		} else {
		verticesList2.push_back(verticesList.at(0));
		verticesList2.push_back(verticesList.at(1));
		verticesList2.push_back(verticesList.at(2));
		verticesList2.push_back(v5p);
		}
		this->buildTetrahedonMesh(verticesList2,newMTetraheds,newMFaces,newMHalfEdges,newMVertices, true);
		verticesList2.clear();

		//TETRAHED 2
		if(i == 1 || i == 2) {
		verticesList2.push_back(verticesList.at(2));
		verticesList2.push_back(v5p);
		verticesList2.push_back(verticesList.at(3));
		verticesList2.push_back(verticesList.at(1));
		} else {
		verticesList2.push_back(verticesList.at(1));
		verticesList2.push_back(v5p);
		verticesList2.push_back(verticesList.at(3));
		verticesList2.push_back(verticesList.at(2));

		}
		this->buildTetrahedonMesh(verticesList2,newMTetraheds,newMFaces,newMHalfEdges,newMVertices, true);
		verticesList2.clear();


		//TETRAHED 3
		if(i == 1 || i == 2) {
		verticesList2.push_back(v5p);
		verticesList2.push_back(verticesList.at(3));
		verticesList2.push_back(verticesList.at(0));
		verticesList2.push_back(verticesList.at(2));

		} else {

		verticesList2.push_back(verticesList.at(0));
		verticesList2.push_back(verticesList.at(3));
		verticesList2.push_back(v5p);
		verticesList2.push_back(verticesList.at(2));
		}

		this->buildTetrahedonMesh(verticesList2,newMTetraheds,newMFaces,newMHalfEdges,newMVertices, true);
		verticesList2.clear();


		//TETRAHED 4
		if(i == 1 || i == 2) {

		verticesList2.push_back(verticesList.at(1));
		verticesList2.push_back(verticesList.at(3));
		verticesList2.push_back(verticesList.at(0));
		verticesList2.push_back(v5p);

		}

		else {
		verticesList2.push_back(verticesList.at(0));
		verticesList2.push_back(verticesList.at(3));
		verticesList2.push_back(verticesList.at(1));
		verticesList2.push_back(v5p);
		}
		this->buildTetrahedonMesh(verticesList2,newMTetraheds,newMFaces,newMHalfEdges,newMVertices, true);
		verticesList2.clear();
		verticesList.clear();

		cout << mTetraheds->size() << endl;

	}

	delete mFaces;
	delete mHalfEdges;
	delete mTetraheds;
	delete mVertices;
	this->mFaces = newMFaces;
	this->mHalfEdges = newMHalfEdges;
	this->mTetraheds = newMTetraheds;
	this->mVertices = newMVertices;


}

arma::Mat<double> TetrahedMesh::pickNode() {


	arma::Mat<double> temp = this->mVertices->at(currentNode).getPosition();

	return temp;
}

void TetrahedMesh::pickNextNode() {

	if (this->currentNode < mVertices->size()-1){
	this->currentNode += 1;

	} else {
	this->currentNode = 0;
	}

}

void TetrahedMesh::pickNextFourNode() {

pickedNodes.clear();

vector<unsigned int > list;
unsigned int currentNode = 0;
arma::Mat<double> temp;
arma::Mat<double> temp2;
temp = mVertices->at(0).getPosition();

    for(int i = 0; i<mVertices->size();i++){

        temp2 = mVertices->at(i).getPosition();
        if(temp(0) > temp2(0)){

        currentNode = i;

        temp = mVertices->at(i).getPosition();

        }

    }

    pickedNodes.push_back(currentNode);

    arma::Mat<double> temp3;
    arma::Mat<double> temp4 = mVertices->at(currentNode).getPosition();
    for(int i = 0; i<mVertices->size();i++){

        temp3 = mVertices->at(i).getPosition();
        if (temp3(0) == temp4(0) && vecEquals(temp3, temp4) == false) {

             pickedNodes.push_back(i);
             if (pickedNodes.size() > 3) {

                break;
             }

        }

    }


}

unsigned int TetrahedMesh::getCurrentNode() {


return this->currentNode;
}


void TetrahedMesh::crackStructure(unsigned int ind, arma::Mat<double> eigVec) {

	Tetrahed tempTetra = mTetraheds->at(ind);

	vector<unsigned int> faceIndices = tempTetra.getFaceInd();
	vector<unsigned int> vertexIndices;

     for(int i=0; i<4;i++){

        Face & face = mFaces->at(faceIndices[i]);

        HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

        vertexIndices.push_back(edge->getVertexInd());

        edge = &mHalfEdges->at(edge->getNextInd());
        vertexIndices.push_back(edge->getVertexInd());

        edge = &mHalfEdges->at(edge->getNextInd());
        vertexIndices.push_back(edge->getVertexInd());
    }

	 int size = vertexIndices.size();

	unsigned int rIndex = rand() % size;
	unsigned int vIndex = vertexIndices.at(rIndex);

	arma::Mat<double> v(3,1);

	v= mVertices->at(vIndex).getPosition();






	vector<unsigned int> adjecentList = getAdjecentTetraheds(ind, v);
	if (adjecentList.size() > 0) {
	determineLocationOfCrack(ind,adjecentList,v,eigVec,vIndex);
	}


}

vector<unsigned int> TetrahedMesh::getAdjecentTetraheds(unsigned int ind, arma::Mat<double> v) {

	vector<unsigned int> tetraIndices;
	Tetrahed tempTetra = mTetraheds->at(ind);
	vector<unsigned int> faceIndices = tempTetra.getFaceInd();
	int count = 0;

	for(int i=0; i<4;i++){

		 //get face
		Face & face = mFaces->at(faceIndices[i]);
		//get all vertices of that face
		HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

		Vertex & v1 = mVertices->at(edge->getVertexInd());
		edge = &mHalfEdges->at(edge->getNextInd());

		Vertex & v2 = mVertices->at(edge->getVertexInd());
		edge = &mHalfEdges->at(edge->getNextInd());

		Vertex & v3 = mVertices->at(edge->getVertexInd());
		int index = face.getOppositeFaceInd();
		if  (index != -1 && (vecEquals(v, v1.getPosition()) || vecEquals(v, v2.getPosition()) || vecEquals(v, v3.getPosition()))){

			getAdjecent(mFaces->at(index).getTetrahedInd(), ind, &tetraIndices, v, count);
			break;
		}


	}
	return tetraIndices;

}

void TetrahedMesh::getAdjecent(unsigned int currentind, unsigned int baseind, vector<unsigned int> *tetraIndices, arma::Mat<double> v, int count){

	count += 1;


	Tetrahed tempTetra = mTetraheds->at(currentind);
	vector<unsigned int> faceIndices = tempTetra.getFaceInd();

	 for(int i=0; i<4;i++){
		bool check = false;
		 //get face
		Face & face = mFaces->at(faceIndices[i]);

		//get all vertices of that face
		HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

		Vertex & v1 = mVertices->at(edge->getVertexInd());
		edge = &mHalfEdges->at(edge->getNextInd());

		Vertex & v2 = mVertices->at(edge->getVertexInd());
		edge = &mHalfEdges->at(edge->getNextInd());

		Vertex & v3 = mVertices->at(edge->getVertexInd());
		int index = face.getOppositeFaceInd();

	if (index != -1 && (vecEquals(v, v1.getPosition()) || vecEquals(v, v2.getPosition()) || vecEquals(v, v3.getPosition()))){


				currentind = mFaces->at(index).getTetrahedInd();
				if(currentind != baseind) {
				 for (int j = 0; j < tetraIndices->size(); j++) {
					 if (tetraIndices->at(j) == currentind) {
						 check = true;
						 break;
					 }
				 }
				 if (check == false) {

				 tetraIndices->push_back(currentind);

				 getAdjecent(currentind,  baseind, tetraIndices, v,count);
				 }


		 }
    }

	}


}


arma::Mat<double> TetrahedMesh::tetCenterOfMass(unsigned int tInd) {

	vector<arma::Mat<double> > vertexList = getVertexPosition(tInd);
	arma::Mat<double> center;
	center = arma::zeros(3,1);
	for(int i = 0; i <vertexList.size(); i++) {

		center += vertexList.at(i).rows(0,2);

	}

	center = center/vertexList.size();

	return center;
}


arma::Mat<double> TetrahedMesh::determineLocationOfCrack(unsigned int tInd, vector<unsigned int> adjecentList, arma::Mat<double> v, arma::Mat<double> eigVec,unsigned int vIndex) {

	double D = arma::dot(trans(v),eigVec);
	vector<unsigned int> Xpos;
	vector<unsigned int> Xneg;

	for(int i = 0; i < adjecentList.size();i++) {

		arma::Mat<double> cMass = tetCenterOfMass(adjecentList.at(i));

		double distanceFromPlane = arma::dot(cMass,eigVec)-D;
		if (distanceFromPlane > 0.0) {

			Xpos.push_back(adjecentList.at(i));

		} else if (distanceFromPlane < 0.0) {

			Xneg.push_back(adjecentList.at(i));
		}
	}

	//deconnect(Xpos, v);
	deconnect(Xneg, v, vIndex);

	arma::Mat<double> mm;
	return mm;
}

void TetrahedMesh::deconnect(vector<unsigned int> adjecentList, arma::Mat<double> v, unsigned int vIndex) {
	unsigned int vInd;
	Vertex temp;
	temp.setPosition(v);

	mVertices->push_back(temp);
	vInd = mVertices->size()-1;

	for (int j = 0; j < adjecentList.size(); j++) {

		vector<unsigned int> faceIndices = mTetraheds->at(adjecentList.at(j)).getFaceInd();

			for (int i = 0; i< faceIndices.size();i++) {

				Face & face = mFaces->at(faceIndices[i]);


				HalfEdge* edge = &mHalfEdges->at(face.getEdgeInd());

				Vertex & v1 = mVertices->at(edge->getVertexInd());

				if(vecEquals(v1.getPosition(),v)) {
					//mFaces->at(faceIndices[i]).setOppositeFaceInd(-1);
					edge->setVertexInd(vInd);

					break;
				}

				edge = &mHalfEdges->at(edge->getNextInd());

				Vertex & v2 = mVertices->at(edge->getVertexInd());

					if(vecEquals(v2.getPosition(),v)) {
					//mFaces->at(faceIndices[i]).setOppositeFaceInd(-1);
					edge->setVertexInd(vInd);


					break;
				}

				edge = &mHalfEdges->at(edge->getNextInd());

				Vertex & v3 = mVertices->at(edge->getVertexInd());

					if(vecEquals(v3.getPosition(),v)) {

					edge->setVertexInd(vInd);


					break;
				}



			}
	}
}


