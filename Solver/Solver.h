#include "armadillo"
#include "ConjugateGradient.h"
#include "../mesh/TetrahedMesh.h"
#include <vector>


class Solver{
public:
    Solver();
    ~Solver();

    void contstructKe(TetrahedMesh *mesh);
    void calcNewPosition(TetrahedMesh *mesh, arma::Mat<double> Fxt);
	void tetrahedronAssemble(arma::Mat<double> &K ,arma::Mat<double> k, int i, int j, int m, int n);
	void setStaticState() {update = true;}
	void checkCollision();

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

};
