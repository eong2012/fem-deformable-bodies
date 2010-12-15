#include "TetrahedMesh.h"
#include "MeshReader.h"
#include <string>

#define INNER   0
#define OUTER   1
#define ALL     2


#ifndef VOLUMEGENERATOR_H
#define VOLUMEGENERATOR_H

class VolumeGenerator {
public:
	VolumeGenerator();
	~VolumeGenerator();
	void generateVolumeFromFile(string filename);
	void generateVolume();
	void createTetra(arma::Mat<double> v1, arma::Mat<double> v2, arma::Mat<double> v3,arma::Mat<double> v4);
	void render();
	void changeNormalRenderMode();
	void changeEdgeRenderMode();
	void changeTriangleRenderMode();

    TetrahedMesh* getTetrahedMesh();
	void subdivide();


private :
	TetrahedMesh *tetrahedMesh;
	MeshReader *meshReader;
	int normalMode;
	int edgeMode;
	int triangleMode;
	vector<arma::Mat<double> > vertices;

};

#endif
