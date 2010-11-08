#include "armadillo"

class ConjugateGradient {
public:
	

	ConjugateGradient(unsigned int maxiter, double tolerance) : mMaxIter(maxiter), mIter(0), mMaxTolerance(tolerance), mTolerance(-1) {}
	unsigned int getNumIter() const { return mIter; }
	unsigned int getMaxNumIter() const { return mMaxIter; }
	double getTolerance() const { return mTolerance; }
	double getMaxTolerance() const { return mMaxTolerance; }

	bool solve(arma::Mat<double> A, arma::Mat<double> *x, arma::Mat<double> b);

	private:
	unsigned int mMaxIter, mIter;
    double mMaxTolerance, mTolerance;
};