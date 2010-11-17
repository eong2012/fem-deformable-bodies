#include "VolumeGenerator.h"


VolumeGenerator::VolumeGenerator() {

	tetrahedMesh = new TetrahedMesh();
	normalMode = 2;
	edgeMode = 2;
	triangleMode = 2;
}

void VolumeGenerator::generateVolume() {

	arma::Mat<double> temp1,temp2,temp3,temp4;
	arma::Mat<double> t = arma::zeros(1,3);
    float scale = 0.7;
	vector<arma::Mat<double> > vertices;

    temp1 << 0.5f << 0.5f << 0.5f;
	temp2 << -0.5f << 0.5f << -0.5f;
    temp3 << 0.5f << -0.5f << -0.5f;
    temp4 << -0.5f << -0.5f << 0.5f;

	vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

    tetrahedMesh->AddTetrahedron(vertices);

	/* Tetrahedrar som GER NEGATIV VOLYM !!11!!11!!!
	vertices.clear();

	temp1 = arma::Mat<double>(0.5f, 0.5f,0.5f);
    temp2 = arma::Mat<double>(-0.5f, -0.5f, 0.5f);
	temp3 = arma::Mat<double>(0.5f, -0.5f, -0.5f);
    temp4 = arma::Mat<double>(0.5f, -0.5f, 0.5f);

    vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

	tetrahedMesh->AddTetrahedron(vertices);
	*//*
	vertices.clear();

	temp1 = arma::Mat<double>(0.5f, 0.5f,0.5f);
    temp2 = arma::Mat<double>(-0.5f, 0.5f,-0.5f);
	temp3 = arma::Mat<double>(0.5f, 0.5f, -0.5f);
	temp4 = arma::Mat<double>(0.5f, -0.5f, -0.5f);


    vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

	tetrahedMesh->AddTetrahedron(vertices);
	*/
	vertices.clear();

	temp1 << 0.5f << 0.5f << 0.5f;
    temp2 << -0.5f << 0.5f <<-0.5f;
	temp3 <<-0.5f << -0.5f << 0.5f;
    temp4 <<-0.5f << 0.5f << 0.5f;

    vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

	tetrahedMesh->AddTetrahedron(vertices);

	vertices.clear();

	temp1 << -0.5f << 0.5f <<-0.5f;
    temp2 << 0.5f << -0.5f << -0.5f;
	temp3 << -0.5f << -0.5f << 0.5f;
    temp4 << -0.5f << -0.5f << -0.5f;

    vertices.push_back(temp1*scale+t);
	vertices.push_back(temp2*scale+t);
	vertices.push_back(temp3*scale+t);
	vertices.push_back(temp4*scale+t);

	tetrahedMesh->AddTetrahedron(vertices);
}

void VolumeGenerator::subdivide() {

	tetrahedMesh->subdivide();

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
