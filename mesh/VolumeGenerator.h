#include "HalfFaceMesh.h"
#include <string>

class VolumeGenerator {
public:
	VolumeGenerator();
	void generateVolumeFromFile(string filename);
	void generateVolume();
	void render();
	

private :
	HalfFaceMesh *halfFaceMesh;

};