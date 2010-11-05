#include "armadillo"
#include "../mesh/TetrahedMesh.h"
#include <vector>


class Solver{
public:
    Solver();
    ~Solver();

    void contstructKe(TetrahedMesh *mesh);
    void calcNewPosition(TetrahedMesh *mesh);

private:

    vector<arma::Mat<double> > mKMatrices;






};
