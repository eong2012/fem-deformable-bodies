#include "TetrahedMesh.h"
#include <string>

#define INNER   0
#define OUTER   1
#define ALL     2


class VolumeGenerator {
public:
	VolumeGenerator();
	void generateVolumeFromFile(string filename);
	void generateVolume();
	void render();
	void changeNormalRenderMode();
	void changeEdgeRenderMode();
	void changeTriangleRenderMode();

	

private :
	TetrahedMesh *tetrahedMesh;
	int normalMode;
	int edgeMode;
	int triangleMode;

};