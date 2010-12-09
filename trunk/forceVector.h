#include "armadillo"

class forceVector {
public: 

	forceVector();
	void translateEndPoint(arma::Mat<double> t);
	arma::Mat<double> startPoint;
	arma::Mat<double> endPoint;


};