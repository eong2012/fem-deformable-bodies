#include "VolumeGenerator.h"


VolumeGenerator::VolumeGenerator() {

	halfFaceMesh = new HalfFaceMesh();
}

void VolumeGenerator::generateVolume() {

	vector<Vector3<float>> vertices;
	
	Vector3<float> temp1(0.5f, 0.5f,0.5f);
	Vector3<float> temp2(-0.5f, 0.5f,-0.5f);
    Vector3<float> temp3(0.5f, -0.5f, -0.5f);
    Vector3<float> temp4(-0.5f, -0.5f, 0.5f);

	vertices.push_back(temp1);
	vertices.push_back(temp2);
	vertices.push_back(temp3);
	vertices.push_back(temp4);

    halfFaceMesh->AddTetrahedron(vertices);
	/*
	vertices.clear();

	Vector3<float> temp12(0.5f, 0.5f,0.5f);
    Vector3<float> temp22(0.5f, -0.5f, -0.5f);
	Vector3<float> temp32(-0.5f, -0.5f, 0.5f);
    Vector3<float> temp42(0.5f, -0.5f, 0.5f);

    vertices.push_back(temp12);
	vertices.push_back(temp22);
	vertices.push_back(temp32);
	vertices.push_back(temp42);

	halfFaceMesh->AddTetrahedron(vertices);
/*
	vertices.clear();

	Vector3<float> temp13(0.5f, 0.5f,0.5f);
    Vector3<float> temp23(0.5f, -0.5f, -0.5f);
	Vector3<float> temp33(-0.5f, 0.5f,-0.5f);
   
    Vector3<float> temp43(0.5f, 0.5f, -0.5f);

    vertices.push_back(temp13);
	vertices.push_back(temp23);
	vertices.push_back(temp33);
	vertices.push_back(temp43);

	halfFaceMesh->AddTetrahedron(vertices);

	vertices.clear();

	Vector3<float> temp14(0.5f, 0.5f,0.5f);
    Vector3<float> temp24(-0.5f, -0.5f, 0.5f);
	Vector3<float> temp34(-0.5f, 0.5f,-0.5f);
    Vector3<float> temp44(-0.5f, 0.5f, 0.5f);

    vertices.push_back(temp14);
	vertices.push_back(temp24);
	vertices.push_back(temp34);
	vertices.push_back(temp44);

	halfFaceMesh->AddTetrahedron(vertices);

	vertices.clear();

	Vector3<float> temp15(-0.5f, 0.5f,-0.5f);
    Vector3<float> temp25(-0.5f, -0.5f, 0.5f);
	Vector3<float> temp35(0.5f, -0.5f, -0.5f);

    Vector3<float> temp45(-0.5f, -0.5f, -0.5f);

    vertices.push_back(temp15);
	vertices.push_back(temp25);
	vertices.push_back(temp35);
	vertices.push_back(temp45);

	halfFaceMesh->AddTetrahedron(vertices);
	
*/
}

void VolumeGenerator::render() {

	halfFaceMesh->Render();
}