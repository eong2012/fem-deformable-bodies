#include "forceVector.h"

forceVector::forceVector() {
	
	this->startPoint = arma::zeros(3,1);
	this->endPoint = arma::zeros(3,1);
	this->endPoint(2) = 50.0;
}

void forceVector::translateEndPoint(arma::Mat<double> t) {

	this->endPoint += t;

}