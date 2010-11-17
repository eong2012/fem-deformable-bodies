#include "VolumeGenerator.h"


VolumeGenerator::VolumeGenerator() {

	tetrahedMesh = new TetrahedMesh();
	normalMode = 2;
	edgeMode = 2;
	triangleMode = 2;
}

void VolumeGenerator::generateVolume() {

	arma::Mat<double> vertex1,vertex2,vertex3,vertex4, vertex5, vertex6, vertex7, vertex8;

	// Skapa ett gäng vertexpunkter

    vertex1 <<  0.5f <<  0.5f <<  0.5f;
	vertex2 << -0.5f <<  0.5f << -0.5f;
    vertex3 <<  0.5f << -0.5f << -0.5f;
    vertex4 << -0.5f << -0.5f <<  0.5f;
    vertex5 <<  0.5f <<  0.5f << -0.5f;
    vertex6 << -0.5f <<  0.5f <<  0.5f;
    vertex7 << -0.5f << -0.5f << -0.5f;
    vertex8 <<  0.5f << -0.5f <<  0.5f;

    //Tetror med negativ volym
    createTetra(vertex1,vertex3, vertex8, vertex4);


    //Tetror med positiv volym
    createTetra(vertex1,vertex2, vertex3, vertex4); //Tetra 1
    createTetra(vertex1,vertex2, vertex4, vertex6); //Tetra 4
    createTetra(vertex2,vertex3, vertex4, vertex7); //Tetra 5
    createTetra(vertex1,vertex2, vertex5, vertex3);


}


void VolumeGenerator::createTetra(arma::Mat<double> v1, arma::Mat<double> v2, arma::Mat<double> v3,arma::Mat<double> v4){
    //Kolla efter positivt volym
    arma::Mat<double> Volmat;
    Volmat = join_rows(arma::ones(4,1),join_cols(join_cols(join_cols(v1,v2),v3),v4));
    double V = det(Volmat)/6.0;
    if(V > 0){
        //Skapa tetras
        arma::Mat<double> t = arma::zeros(1,3);
        float scale = 1.0;

        vertices.clear();
        vertices.push_back(v1*scale+t);
        vertices.push_back(v2*scale+t);
        vertices.push_back(v3*scale+t);
        vertices.push_back(v4*scale+t);
        std::cout << V << endl;

        tetrahedMesh->AddTetrahedron(vertices);
    }
    else
        std::cout << V << endl;


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

