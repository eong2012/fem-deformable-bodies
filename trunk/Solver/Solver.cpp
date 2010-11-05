#include "Solver.h"


Solver::Solver(){

}

Solver::~Solver(){}

void Solver::contstructKe(TetrahedMesh *mesh){



    //for each vertex in tetraheder, get positon
    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();
    vector<Vector3<float> > vPositions;
    for(int k = 0; k < nrOfTetraheds; k++) {

            vPositions = mesh->getVertexPosition(k);

            arma::Col<double> x0(3);
            arma::Col<double> x1(3);
            arma::Col<double> x2(3);
            arma::Col<double> x3(3);

            x0(0) = vPositions[0][0];
            x0(1) = vPositions[0][1];
            x0(2) = vPositions[0][2];

            x1(0) = vPositions[1][0];
            x1(1) = vPositions[1][1];
            x1(2) = vPositions[1][2];

            x2(0) = vPositions[2][0];
            x2(1) = vPositions[2][1];
            x2(2) = vPositions[2][2];

            x3(0) = vPositions[3][0];
            x3(1) = vPositions[3][1];
            x3(2) = vPositions[3][2];


            arma::Mat<double> X;

            X = join_rows(join_rows((x1-x0),(x2-x0)), (x3-x0));

            double V = det(X);

            X = inv(X);

            //beräknar motorseglare :D
            arma::Row<double> y1 = X.row(0);
            arma::Row<double> y2 = X.row(1);
            arma::Row<double> y3 = X.row(2);
            arma::Row<double> y0 = -y1 -y2 -y3;

            vector<arma::Row<double> > y;


            y.push_back(y0);
            y.push_back(y1);
            y.push_back(y2);
            y.push_back(y3);

            float a, b, c, vn, E;
            vn = 0.3;
            E = 0.05;
            //Paranthesis overflow :X DONT DIVIDE BY ZERO
            a = V*E*((1-vn)/((1+vn)*(1-2*vn)));
            b = V*E*(vn/((1+vn)*(1-2*vn)));
            c = V*E*((1-2*vn)/((1+vn)*(1-2*vn)));

            arma::Mat<double> AB;
            arma::Mat<double> C;
            AB << a << b << b << arma::endr
               << b << a << b << arma::endr
               << b << b << a << arma::endr;

            C  << c << 0 << 0 << arma::endr
               << 0 << c << 0 << arma::endr
               << 0 << 0 << c << arma::endr;

            arma::Mat<double> K1,K2,K3,K4, Ksub;
            arma::Mat<double> Ke(12,12);
            Ke = arma::zeros<arma::mat>(12,12);

            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){

                    K1 << y[i][0] << 0 << 0 << arma::endr
                          << 0 << y[i][1] << 0 << arma::endr
                          << 0 << 0 << y[i][2] << arma::endr;
                    K2 << y[j][0] << 0 << 0 << arma::endr
                          << 0 << y[j][1] << 0 << arma::endr
                          << 0 << 0 << y[j][2] << arma::endr;
                    K3 << y[i][1] << 0 << y[i][2] << arma::endr
                         << y[i][0] << y[i][2] << 0 << arma::endr
                          << 0 << y[i][1] << y[i][1] << arma::endr;
                    K4 << y[j][1] << 0 << y[j][2] << arma::endr
                          << y[j][0] << y[j][2] << 0 << arma::endr
                          << 0 << y[j][1] << y[j][1] << arma::endr;

                    K4 = trans(K4); //Same as Ksub3 but transp

                    //.submat( first_row, first_col, last_row, last_col )
                    Ksub = K1*AB*K2 + K3*C*K4;

                    Ke.submat(i*3,j*3,i*3+2,j*3+2) =  K1*AB*K2 + K3*C*K4;
                }
            }

            this->mKMatrices.push_back(Ke);
            vPositions.clear();

    }


    //construct  y

    //calculate a, b, c

    //define matrices

    //k_e = matrix_yi * matrix_ab * matrix_yj +  matrix_yi2 * matrix_c * matrix_yi2


}

void Solver::calcNewPosition(TetrahedMesh *mesh)
{
     //for each vertex in tetraheder, get positon
    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();
    vector<Vector3<float> > vPositions;



    arma::mat u = arma::randu<arma::mat>(12,1)*100.0;


    for(int k = 0; k < nrOfTetraheds; k++) {

            vPositions = mesh->getVertexPosition(k);

            arma::Col<double> x0(3);
            arma::Col<double> x1(3);
            arma::Col<double> x2(3);
            arma::Col<double> x3(3);

            x0(0) = vPositions[0][0];
            x0(1) = vPositions[0][1];
            x0(2) = vPositions[0][2];

            x1(0) = vPositions[1][0];
            x1(1) = vPositions[1][1];
            x1(2) = vPositions[1][2];

            x2(0) = vPositions[2][0];
            x2(1) = vPositions[2][1];
            x2(2) = vPositions[2][2];

            x3(0) = vPositions[3][0];
            x3(1) = vPositions[3][1];
            x3(2) = vPositions[3][2];

            //Store the vertices in a 12x1 vector
            arma::Mat<double> x = join_cols(join_cols(join_cols(x0,x1), x2),x3);

            arma::Mat<double> force = mKMatrices[k]*u;

            double dt = 0.01;

            //euler

				arma::Mat<double> v(12,1);
				arma::Mat<double> forceExt(12,1);
				forceExt = arma::ones(12,1)*0;
				
            
				v = dt*(force+forceExt);
			
			
		
				
				
			

            arma::Mat<double> newX(12,1);
            newX = x + dt*v;

            vector<Vector3<float> > tempX;
            for(int i = 0; i < 4;i++){
                Vector3<float> tx(newX(i*3),newX(i*3+1), newX(i*3+2));
                tempX.push_back(tx);
            }

            mesh->updatePosition(k,tempX);


            //clear positions vector for the new tetra
            vPositions.clear();


    }
}
