#include "Solver.h"


Solver::Solver(){

	conjugateGradient = new ConjugateGradient(10.0,0.00000000002);
	K = arma::zeros(25,25);
	Xpre = arma::zeros(24,1);
	Vpre = arma::zeros(24,1);
	check = false;
	update = false;
	grav = arma::zeros(24,1);
	int i = 0;
	while(i*3+1 < grav.n_rows) {
	
	grav(i*3+1) = -9.82;
	i++;
	}
}

Solver::~Solver(){

	delete[] this->vOVer1;
	delete[] this->vOVer2;
	delete[] this->xOVer1;
	delete[] this->xOVer2;
}

void Solver::contstructKe(TetrahedMesh *mesh){



    //for each vertex in tetraheder, get positon
    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();

	
	xOVer1 = new arma::Mat<double>[nrOfTetraheds];
	xOVer2 = new arma::Mat<double>[nrOfTetraheds];
	vOVer1 = new arma::Mat<double>[nrOfTetraheds];
	vOVer2 = new arma::Mat<double>[nrOfTetraheds];

		for(int i = 0; i < nrOfTetraheds; i++){
		
			xOVer1[i] = arma::zeros(12,1);
			xOVer2[i] = arma::zeros(12,1);
			vOVer1[i] = arma::zeros(12,1);
			vOVer2[i] = arma::zeros(12,1);
		
		}

    vector<arma::Mat<double> > vPositions;
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

			arma::Mat<double> xOrgin = join_cols(join_cols(join_cols(x0,x1), x2),x3);
			this->xOrgin.push_back(xOrgin);
	
            arma::Mat<double> X, Volmat;

            X = join_rows(join_rows((x1-x0),(x2-x0)), (x3-x0));
			
			
			Volmat = join_rows(arma::ones(4,1),join_cols(join_cols(join_cols(trans(x0),trans(x1)),trans(x2)),trans(x3)));
			
            double V = det(Volmat)/6.0;
			

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
            vn = 0.6;
            E = 0.8;
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
			cout << vPositions[0][3] << " " << vPositions[1][3] << " " << vPositions[2][3] << " " << vPositions[3][3] << endl;
			this->tetrahedronAssemble(this->K,Ke,vPositions[0][3]+1,vPositions[1][3]+1,vPositions[2][3]+1, vPositions[3][3]+1);
            this->mKMatrices.push_back(Ke);

            vPositions.clear();

    }


    //construct  y

    //calculate a, b, c

    //define matrices

    //k_e = matrix_yi * matrix_ab * matrix_yj +  matrix_yi2 * matrix_c * matrix_yi2
		cout << "------" << endl;
		K.shed_col(0);
		K.shed_row(0);
	//cout << this->K << endl;

}

void Solver::calcNewPosition(TetrahedMesh *mesh, arma::Mat<double> Fxt)
{
     //for each vertex in tetraheder, get positon
    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();
    vector<arma::Mat<double> > vPositions;
	arma::Mat<double> M = arma::zeros<arma::mat>(24,24);
	M = arma::eye(24,24)*1;

	arma::Mat<double> X(24,1);
	X(0) = mesh->mVertices[0].getPosition()[0];
	X(1) = mesh->mVertices[0].getPosition()[1];
	X(2) = mesh->mVertices[0].getPosition()[2];
	X(3) = mesh->mVertices[1].getPosition()[0];
	X(4) = mesh->mVertices[1].getPosition()[1];
	X(5) = mesh->mVertices[1].getPosition()[2];
	X(6) = mesh->mVertices[2].getPosition()[0];
	X(7) = mesh->mVertices[2].getPosition()[1];
	X(8) = mesh->mVertices[2].getPosition()[2];
	X(9) = mesh->mVertices[3].getPosition()[0];
	X(10) = mesh->mVertices[3].getPosition()[1];
	X(11) = mesh->mVertices[3].getPosition()[2];
	X(12) = mesh->mVertices[4].getPosition()[0];
	X(13) = mesh->mVertices[4].getPosition()[1];
	X(14) = mesh->mVertices[4].getPosition()[2];
	X(15) = mesh->mVertices[5].getPosition()[0];
	X(16) = mesh->mVertices[5].getPosition()[1];
	X(17) = mesh->mVertices[5].getPosition()[2];
	X(18) = mesh->mVertices[6].getPosition()[0];
	X(19) = mesh->mVertices[6].getPosition()[1];
	X(20) = mesh->mVertices[6].getPosition()[2];
	X(21) = mesh->mVertices[7].getPosition()[0];
	X(22) = mesh->mVertices[7].getPosition()[1];
	X(23) = mesh->mVertices[7].getPosition()[2];


    //arma::mat u = arma::randu<arma::mat>(24,1)*0.0;
	//u(0) = -0.0001;
	//u(23) = -0.0001;
	//arma::mat fex = arma::zeros(24,1);
	//fex(0) = 0.001;
	//fex(1) = 0.001;
	//fex(2) = 0.001;

	if (check == false) {
		Xpre = X;
		check = true;
	} 
	
	if(update == true) {
		Xpre = X;
		
	}

	arma::Mat<double> u;
	u = X-Xpre;
	
	//arma::norm(u,1);

	float dt = 0.01;
	arma::Mat<double> force;
	force = this->K*(u)+Fxt+grav*0.0;
	cout << u << endl;
	
	//cout << force << endl;
	//arma::Mat<double> b = arma::ones(24,1);
	
	//b = (M*(X-X) + this->K*(u));
	
	
	//cout << X-Xpre << endl;
	arma::Mat<double> v;
	v = Vpre+dt*force;
	v = v;
	Vpre = v;
	
	arma::Mat<double> x = X+dt*v;
	

	

	mesh->mVertices[0].setPosition(Vector3<float>(x(0),x(1), x(2)));
	mesh->mVertices[1].setPosition(Vector3<float>(x(3),x(4), x(5)));
	mesh->mVertices[2].setPosition(Vector3<float>(x(6),x(7), x(8)));
	mesh->mVertices[3].setPosition(Vector3<float>(x(9),x(10), x(11)));
	mesh->mVertices[4].setPosition(Vector3<float>(x(12),x(13), x(14)));
	mesh->mVertices[5].setPosition(Vector3<float>(x(15),x(16), x(17)));
	mesh->mVertices[6].setPosition(Vector3<float>(x(18),x(19), x(20)));
	mesh->mVertices[7].setPosition(Vector3<float>(x(21),x(22), x(23)));

	//cout << nrOfTetraheds << endl;
/*
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


			
			v = dt*f;

			/*
            //Store the vertices in a 12x1 vector
            arma::Mat<double> x = join_cols(join_cols(join_cols(x0,x1), x2),x3);
		
			if(k % 2 == 0) {
			xOVer2[k] = x;
			} else {
			xOVer1[k] = x ;
			}
            //arma::Mat<double> force = mKMatrices[k]*u;
		
            double dt = 0.1;
			
			arma::Mat<double> *v;
			v = new arma::Mat<double>(12,1);

			arma::Mat<double> b(12,1);
			
			b = arma::ones(12,1);
			if(k % 2 == 0) {
				b = M*vOVer1[k]+dt*(mKMatrices[k]*(x-xOVer2[k]));
			} 
			else {
				b = M*vOVer2[k]+dt*(mKMatrices[k]*(x-xOVer1[k]));
			}
            //euler
			bool solved = conjugateGradient->solve(M-dt*dt*mKMatrices[k], v, b);
			//cout << solved << endl;
			
			
			if(k % 2 == 0) {
			vOVer2[k] = (*v);
			} else {
			vOVer1[k] = (*v);
			
			}
			
			
		    arma::Mat<double> newX(12,1);
            if(k % 2 == 0) {
				newX = x+(vOVer2[k]*dt);
				
			} else {

			newX = x+(vOVer1[k]*dt);
			
			
			}
			
			
			
            vector<Vector3<float> > tempX;
            for(int i = 0; i < 4;i++){
                Vector3<float> tx(newX(i*3),newX(i*3+1), newX(i*3+2));
                tempX.push_back(tx);
            }

            mesh->updatePosition(k,tempX);
			

            //clear positions vector for the new tetra
            vPositions.clear();
		*/

    //}
}

void Solver::tetrahedronAssemble(arma::Mat<double> &K ,arma::Mat<double> k, int i, int j, int m, int n){
	
	int offset = -1;
	K(3*i-2,3*i-2) = K(3*i-2,3*i-2) + k(1+offset,1+offset);
	K(3*i-2,3*i-1) = K(3*i-2,3*i-1) + k(1+offset,2+offset);
	K(3*i-2,3*i) = K(3*i-2,3*i) + k(1+offset,3+offset);
	K(3*i-2,3*j-2) = K(3*i-2,3*j-2) + k(1+offset,4+offset);
	K(3*i-2,3*j-1) = K(3*i-2,3*j-1) + k(1+offset,5+offset);
	K(3*i-2,3*j) = K(3*i-2,3*j) + k(1+offset,6+offset);
	K(3*i-2,3*m-2) = K(3*i-2,3*m-2) + k(1+offset,7+offset);
	K(3*i-2,3*m-1) = K(3*i-2,3*m-1) + k(1+offset,8+offset);
	K(3*i-2,3*m) = K(3*i-2,3*m) + k(1+offset,9+offset);
	K(3*i-2,3*n-2) = K(3*i-2,3*n-2) + k(1+offset,10+offset);
	K(3*i-2,3*n-1) = K(3*i-2,3*n-1) + k(1+offset,11+offset);
	K(3*i-2,3*n) = K(3*i-2,3*n) + k(1+offset,12+offset);
	K(3*i-1,3*i-2) = K(3*i-1,3*i-2) + k(2+offset,1+offset);
	K(3*i-1,3*i-1) = K(3*i-1,3*i-1) + k(2+offset,2+offset);
	K(3*i-1,3*i) = K(3*i-1,3*i) + k(2+offset,3+offset);
	K(3*i-1,3*j-2) = K(3*i-1,3*j-2) + k(2+offset,4+offset);
	K(3*i-1,3*j-1) = K(3*i-1,3*j-1) + k(2+offset,5+offset);
	K(3*i-1,3*j) = K(3*i-1,3*j) + k(2+offset,6+offset);
	K(3*i-1,3*m-2) = K(3*i-1,3*m-2) + k(2+offset,7+offset);
	K(3*i-1,3*m-1) = K(3*i-1,3*m-1) + k(2+offset,8+offset);
	K(3*i-1,3*m) = K(3*i-1,3*m) + k(2+offset,9+offset);
	K(3*i-1,3*n-2) = K(3*i-1,3*n-2) + k(2+offset,10+offset);
	K(3*i-1,3*n-1) = K(3*i-1,3*n-1) + k(2+offset,11+offset);
	K(3*i-1,3*n) = K(3*i-1,3*n) + k(2+offset,12+offset);
	K(3*i,3*i-2) = K(3*i,3*i-2) + k(3+offset,1+offset);
	K(3*i,3*i-1) = K(3*i,3*i-1) + k(3+offset,2+offset);
	K(3*i,3*i) = K(3*i,3*i) + k(3+offset,3+offset);
	K(3*i,3*j-2) = K(3*i,3*j-2) + k(3+offset,4+offset);
	K(3*i,3*j-1) = K(3*i,3*j-1) + k(3+offset,5+offset);
	K(3*i,3*j) = K(3*i,3*j) + k(3+offset,6+offset);
	K(3*i,3*m-2) = K(3*i,3*m-2) + k(3+offset,7+offset);
	K(3*i,3*m-1) = K(3*i,3*m-1) + k(3+offset,8+offset);
	K(3*i,3*m) = K(3*i,3*m) + k(3+offset,9+offset);
	K(3*i,3*n-2) = K(3*i,3*n-2) + k(3+offset,10+offset);
	K(3*i,3*n-1) = K(3*i,3*n-1) + k(3+offset,11+offset);
	K(3*i,3*n) = K(3*i,3*n) + k(3+offset,12+offset);
	K(3*j-2,3*i-2) = K(3*j-2,3*i-2) + k(4+offset,1+offset);
	K(3*j-2,3*i-1) = K(3*j-2,3*i-1) + k(4+offset,2+offset);
	K(3*j-2,3*i) = K(3*j-2,3*i) + k(4+offset,3+offset);
	K(3*j-2,3*j-2) = K(3*j-2,3*j-2) + k(4+offset,4+offset);
	K(3*j-2,3*j-1) = K(3*j-2,3*j-1) + k(4+offset,5+offset);
	K(3*j-2,3*j) = K(3*j-2,3*j) + k(4+offset,6+offset);
	K(3*j-2,3*m-2) = K(3*j-2,3*m-2) + k(4+offset,7+offset);
	K(3*j-2,3*m-1) = K(3*j-2,3*m-1) + k(4+offset,8+offset);
	K(3*j-2,3*m) = K(3*j-2,3*m) + k(4+offset,9+offset);
	K(3*j-2,3*n-2) = K(3*j-2,3*n-2) + k(4+offset,10+offset);
	K(3*j-2,3*n-1) = K(3*j-2,3*n-1) + k(4+offset,11+offset);
	K(3*j-2,3*n) = K(3*j-2,3*n) + k(4+offset,12+offset);
	K(3*j-1,3*i-2) = K(3*j-1,3*i-2) + k(5+offset,1+offset);
	K(3*j-1,3*i-1) = K(3*j-1,3*i-1) + k(5+offset,2+offset);
	K(3*j-1,3*i) = K(3*j-1,3*i) + k(5+offset,3+offset);
	K(3*j-1,3*j-2) = K(3*j-1,3*j-2) + k(5+offset,4+offset);
	K(3*j-1,3*j-1) = K(3*j-1,3*j-1) + k(5+offset,5+offset);
	K(3*j-1,3*j) = K(3*j-1,3*j) + k(5+offset,6+offset);
	K(3*j-1,3*m-2) = K(3*j-1,3*m-2) + k(5+offset,7+offset);
	K(3*j-1,3*m-1) = K(3*j-1,3*m-1) + k(5+offset,8+offset);
	K(3*j-1,3*m) = K(3*j-1,3*m) + k(5+offset,9+offset);
	K(3*j-1,3*n-2) = K(3*j-1,3*n-2) + k(5+offset,10+offset);
	K(3*j-1,3*n-1) = K(3*j-1,3*n-1) + k(5+offset,11+offset);
	K(3*j-1,3*n) = K(3*j-1,3*n) + k(5+offset,12+offset);
	K(3*j,3*i-2) = K(3*j,3*i-2) + k(6+offset,1+offset);
	K(3*j,3*i-1) = K(3*j,3*i-1) + k(6+offset,2+offset);
	K(3*j,3*i) = K(3*j,3*i) + k(6+offset,3+offset);
	K(3*j,3*j-2) = K(3*j,3*j-2) + k(6+offset,4+offset);
	K(3*j,3*j-1) = K(3*j,3*j-1) + k(6+offset,5+offset);
	K(3*j,3*j) = K(3*j,3*j) + k(6+offset,6+offset);
	K(3*j,3*m-2) = K(3*j,3*m-2) + k(6+offset,7+offset);
	K(3*j,3*m-1) = K(3*j,3*m-1) + k(6+offset,8+offset);
	K(3*j,3*m) = K(3*j,3*m) + k(6+offset,9+offset);
	K(3*j,3*n-2) = K(3*j,3*n-2) + k(6+offset,10+offset);
	K(3*j,3*n-1) = K(3*j,3*n-1) + k(6+offset,11+offset);
	K(3*j,3*n) = K(3*j,3*n) + k(6+offset,12+offset);
	K(3*m-2,3*i-2) = K(3*m-2,3*i-2) + k(7+offset,1+offset);
	K(3*m-2,3*i-1) = K(3*m-2,3*i-1) + k(7+offset,2+offset);
	K(3*m-2,3*i) = K(3*m-2,3*i) + k(7+offset,3+offset);
	K(3*m-2,3*j-2) = K(3*m-2,3*j-2) + k(7+offset,4+offset);
	K(3*m-2,3*j-1) = K(3*m-2,3*j-1) + k(7+offset,5+offset);
	K(3*m-2,3*j) = K(3*m-2,3*j) + k(7+offset,6+offset);
	K(3*m-2,3*m-2) = K(3*m-2,3*m-2) + k(7+offset,7+offset);
	K(3*m-2,3*m-1) = K(3*m-2,3*m-1) + k(7+offset,8+offset);
	K(3*m-2,3*m) = K(3*m-2,3*m) + k(7+offset,9+offset);
	K(3*m-2,3*n-2) = K(3*m-2,3*n-2) + k(7+offset,10+offset);
	K(3*m-2,3*n-1) = K(3*m-2,3*n-1) + k(7+offset,11+offset);
	K(3*m-2,3*n) = K(3*m-2,3*n) + k(7+offset,12+offset);
	K(3*m-1,3*i-2) = K(3*m-1,3*i-2) + k(8+offset,1+offset);
	K(3*m-1,3*i-1) = K(3*m-1,3*i-1) + k(8+offset,2+offset);
	K(3*n,3*n-1) = K(3*n,3*n-1) + k(12+offset,11+offset);
	K(3*n,3*n) = K(3*n,3*n) + k(12+offset,12+offset);


}
