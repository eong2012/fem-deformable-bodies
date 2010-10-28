#include "HalfFaceMesh.h"
HalfFaceMesh::HalfFaceMesh(){}
HalfFaceMesh::~HalfFaceMesh(){}

void HalfFaceMesh::setTetraGeometry()
{
    Vector3<float> temp1(0.0f,0.0f,0.0f);
    Vector3<float> temp2(1.0f,0.0f,0.0f);
    Vector3<float> temp3(0.5f,1.0f,5.0f);
    Vector3<float> temp4(1.0f,0.0f,1.0f);

    unsigned int index1,index2,index3,index4;
    bool success = true;
    success &= AddVertex(temp1,index1);
    success &= AddVertex(temp2,index2);
    success &= AddVertex(temp3,index3);
    success &= AddVertex(temp4,index4);


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
