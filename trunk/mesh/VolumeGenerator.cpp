#include "VolumeGenerator.h"


VolumeGenerator::VolumeGenerator() {

	tetrahedMesh = new TetrahedMesh();
	normalMode = 2;
	edgeMode = 2;
	triangleMode = 2;
}

void VolumeGenerator::generateVolume() {
	
	Vector3<float> temp1,temp2,temp3,temp4;
	Vector3<float> t =  Vector3<float>(0.0f, 0.0f,0.0f);
    float scale = 0.7;
	vector<Vector3<float> > vertices;

	temp1 = Vector3<float>(0.5f, 0.5f,0.5f);
	temp2 = Vector3<float>(-0.5f, 0.5f,-0.5f);
    temp3 = Vector3<float>(0.5f, -0.5f, -0.5f);
    temp4 = Vector3<float>(-0.5f, -0.5f, 0.5f);

	vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

    tetrahedMesh->AddTetrahedron(vertices);

	vertices.clear();

	temp1 = Vector3<float>(0.5f, 0.5f,0.5f);
    temp2 = Vector3<float>(-0.5f, -0.5f, 0.5f);
	temp3 = Vector3<float>(0.5f, -0.5f, -0.5f);
    temp4 = Vector3<float>(0.5f, -0.5f, 0.5f);

    vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

	tetrahedMesh->AddTetrahedron(vertices);

	vertices.clear();

	temp1 = Vector3<float>(0.5f, 0.5f,0.5f);
    temp2 = Vector3<float>(-0.5f, 0.5f,-0.5f);
	temp3 = Vector3<float>(0.5f, 0.5f, -0.5f);
	temp4 = Vector3<float>(0.5f, -0.5f, -0.5f);


    vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

	tetrahedMesh->AddTetrahedron(vertices);

	vertices.clear();

	temp1 = Vector3<float>(0.5f, 0.5f,0.5f);
    temp2 = Vector3<float>(-0.5f, 0.5f,-0.5f);
	temp3 = Vector3<float>(-0.5f, -0.5f, 0.5f);
    temp4 = Vector3<float>(-0.5f, 0.5f, 0.5f);

    vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

	tetrahedMesh->AddTetrahedron(vertices);

	vertices.clear();

	temp1 = Vector3<float>(-0.5f, 0.5f,-0.5f);
    temp2 = Vector3<float>(0.5f, -0.5f, -0.5f);
	temp3= Vector3<float>(-0.5f, -0.5f, 0.5f);
    temp4 = Vector3<float>(-0.5f, -0.5f, -0.5f);

    vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

	tetrahedMesh->AddTetrahedron(vertices);


}


void VolumeGenerator::changeNormalRenderMode() {
	if (normalMode > 2 || normalMode < 0 )
	{
		normalMode = 2;

	}else {
	normalMode--;
	}
}

void VolumeGenerator::changeEdgeRenderMode() {
	if (edgeMode > 2 || edgeMode < 0 )
	{
		edgeMode = 2;
	} else {
	edgeMode--;
	}
}

void VolumeGenerator::changeTriangleRenderMode(){
	if (triangleMode > 2 || triangleMode < 0 )
	{
		triangleMode = 2;
	} else {
	triangleMode--;
	}
}


void changeEdgeRenderMode(){

}

void VolumeGenerator::render() {

	tetrahedMesh->Render(triangleMode);
	tetrahedMesh->RenderNormals(normalMode);
	tetrahedMesh->RenderEdges(edgeMode);
}

TetrahedMesh* VolumeGenerator::getTetrahedMesh()
{
    return tetrahedMesh;
}
