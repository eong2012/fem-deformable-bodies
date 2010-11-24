
#include "ConjugateGradient.h"



bool ConjugateGradient::solve(arma::Mat<double> A, arma::Mat<double> *x, arma::Mat<double> b){
    double resid;
    arma::Mat<double> p, q;
    double alpha, beta, rho, rho_1;
    unsigned int max_iter = mMaxIter;
    double tol = mMaxTolerance;

	double normb = arma::norm(b, 2);

    arma::Mat<double> r = b - A*(*x); // N + M*N*fill

    if (normb == 0.0)
      normb = 1;
    if ((resid = arma::norm(r, 2) / normb) <= tol) { // N
      mTolerance = resid;
      mIter = 0;
      return true;
    }

    for (unsigned int i = 1; i <= max_iter; i++) {
      rho = dot(r, r); // N
      if (i == 1)
        p = r; // N
      else {
        beta = rho / rho_1;
        p = r + beta * p; // N + N + N
      }

      q = A*p; // N + M*N*fill
      alpha = rho / dot(p, q); // N

      (*x) = (*x) + alpha * p; // N + N
      r = r - alpha * q; // N + N

      if ((resid = arma::norm(b, 2) / normb) <= tol) { // N
        mTolerance = resid;
        mIter = i;

        return true;
      }
      rho_1 = rho;
    }

    mTolerance = resid;
    mIter = mMaxIter;
    return false;
  }
