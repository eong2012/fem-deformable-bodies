#include "../../usr/include/armadillo.h"
#include "ConjugateGradient.h"
#include "../mesh/TetrahedMesh.h"
#include <vector>


class Solver{
public:
	Solver();
	Solver(int nrOfNodes);
    ~Solver();

    void contstructKe(TetrahedMesh *mesh);
    void calcNewPosition(TetrahedMesh *mesh, arma::Mat<double> Fxt);
	void tetrahedronAssemble(arma::Mat<double> &K ,arma::Mat<double> k, int i, int j, int m, int n);
	void setStaticState() {update = true;}
	void planeCollisionDetection(arma::Mat<double> X);
	void planeCollisionHandler(unsigned int forceIndex);

private:

	ConjugateGradient *conjugateGradient;
    vector<arma::Mat<double> > mKMatrices;
	arma::Mat<double> K;
	arma::Mat<double>  *xOVer1;
	arma::Mat<double>  *xOVer2;
	arma::Mat<double>  *vOVer1;
	arma::Mat<double>  *vOVer2;
	bool check,update;

	arma::Mat<double> Vpre;
	arma::Mat<double> Xpre;
	arma::Mat<double> grav;
	vector<arma::Mat<double> >  xOrgin;
	arma::Mat<double> normal;
	int nrOfNodes;
	arma::Mat<double> ForcePrev;

	arma::Mat<double> v;
	arma::Mat<double> X;
	arma::Mat<double> collisionForce;
	double dt;
	double mass;

	vector<Vertex>* mOriginalPos;
	arma::Mat<double> localVpre;

};