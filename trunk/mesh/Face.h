

#include "../../usr/include/armadillo.h"

#ifndef FACE_H
#define FACE_H

class Face{
public:
    Face();
    ~Face();

    void setNormal(arma::Mat<double> n);
    void setEdgeInd(unsigned int e);
	void setOppositeFaceInd(int f);
    arma::Mat<double> getNormal();
    unsigned int getEdgeInd();
	int getOppositeFaceInd();


private:
    arma::Mat<double> normal;
    unsigned int edgeInd;
	int oppositeFaceInd;

};

#endif
