#include "../../usr/include/armadillo.h"
#include "ConjugateGradient.h"
#include "../mesh/TetrahedMesh.h"
#include <vector>


class Solver{
public:
	Solver();
	Solver(int nrOfNodes);
    ~Solver();

    void constructKe(TetrahedMesh *mesh);
    void calcNewPosition(TetrahedMesh *mesh, arma::Mat<double> Fxt);
	void tetrahedronAssemble(arma::Mat<double> &K ,arma::Mat<double> k, int i, int j, int m, int n);
	arma::Mat<double> calculateB(arma::Mat<double> x1,arma::Mat<double> x2, arma::Mat<double> x3, arma::Mat<double> x4);
	arma::Mat<double> calculateD(float E,float NU);

	arma::Mat<double> calculateRotation(arma::Mat<double> X, arma::Mat<double> X1);
	arma::Mat<double> findRotation(TetrahedMesh *mesh,unsigned int theInd);
	arma::Mat<double> TetrahedronElementStiffness(float E,float NU,arma::Mat<double> x1,arma::Mat<double> x2, arma::Mat<double> x3, arma::Mat<double> x4);
	arma::Mat<double> calculateStress(TetrahedMesh *mesh,unsigned int theInd, float E,float NU, arma::Mat<double> R);
	arma::Mat<double> calculateLargestEIG(arma::Mat<double> stresstensor);
	arma::Mat<double> crackIT(TetrahedMesh *mesh, unsigned int theInd, arma::Mat<double> stresstensor);
	void changeFracture();

	void planeCollisionDetection(arma::Mat<double> X);
	void planeCollisionHandler(unsigned int forceIndex);
	void constructMe(TetrahedMesh *mesh);

	void setParameter(float mass, float fracture);


private:

	ConjugateGradient *conjugateGradient;
    vector<arma::Mat<double> > mKMatrices;
	vector<arma::Mat<double> > mBMatrices;

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

	arma::Mat<double> K;
	arma::Mat<double> C;
	arma::Mat<double> M;

	arma::Mat<double> stress;

	double dt;
	double mass;
	double density;

	double E;
	double alpha, beta;
	double vn;

	bool allowFracture;
	double FractureThresh;

	vector<Vertex>* mOriginalPos;
	arma::Mat<double> localVpre;

};
