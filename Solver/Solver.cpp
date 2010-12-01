#include "Solver.h"


Solver::Solver(int nrOfNodes){

	this->nrOfNodes = nrOfNodes;
	conjugateGradient = new ConjugateGradient(100,0.001);
	K = arma::zeros(this->nrOfNodes*3+1,this->nrOfNodes*3+1);
	M = arma::zeros(this->nrOfNodes*3+1,this->nrOfNodes*3+1);
	Xpre = arma::zeros(this->nrOfNodes*3,1);
	Vpre = arma::zeros(this->nrOfNodes*3,1);
	localVpre = arma::zeros(this->nrOfNodes*3,1);
	check = false;
	update = false;
	grav = arma::zeros(this->nrOfNodes*3,1);
	unsigned int i = 0;
	X = arma::zeros(this->nrOfNodes*3,1);
	dt = 0.001;

	density =1.0;
	this->v = arma::zeros(this->nrOfNodes*3,1);

	collisionForce = arma::zeros(this->nrOfNodes*3,1);

	while(i*3+2 < grav.n_rows) {

	grav(i*3+2) = 9.82*0;
	i++;
	}

	ForcePrev = arma::zeros(this->nrOfNodes*3,1);


}

Solver::~Solver(){

	delete[] this->vOVer1;
	delete[] this->vOVer2;
	delete[] this->xOVer1;
	delete[] this->xOVer2;
}

void Solver::constructKe(TetrahedMesh *mesh){

    mOriginalPos = new vector<Vertex>((*mesh->mVertices));

    //for each vertex in tetraheder, get positon
    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();


	xOVer1 = new arma::Mat<double>[nrOfTetraheds];
	xOVer2 = new arma::Mat<double>[nrOfTetraheds];
	vOVer1 = new arma::Mat<double>[nrOfTetraheds];
	vOVer2 = new arma::Mat<double>[nrOfTetraheds];

		for(unsigned int i = 0; i < nrOfTetraheds; i++){

			xOVer1[i] = arma::zeros(12,1);
			xOVer2[i] = arma::zeros(12,1);
			vOVer1[i] = arma::zeros(12,1);
			vOVer2[i] = arma::zeros(12,1);

		}

    vector<arma::Mat<double> > vPositions;

    for(unsigned int k = 0; k < nrOfTetraheds; k++) {

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
			cout << endl << "Volume: " << V << endl;
            mass += V*density;

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
            vn = 0.33;
            E = 5000.300000;

            //Paranthesis overflow :X DONT DIVIDE BY ZERO
			/*
            a = V*E*((1-vn)/((1+vn)*(1-2*vn)));
            b = V*E*(vn/((1+vn)*(1-2*vn)));
            c = V*E*((1-2*vn)/(2*((1+vn)*(1-2*vn))));

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
            }*/
			arma::Mat<double> Ke(12,12);
			Ke = arma::zeros<arma::mat>(12,12);
			Ke = TetrahedronElementStiffness(E,vn, x0, x1, x2, x3);
			
			this->tetrahedronAssemble(this->K,Ke,vPositions[0][3]+1,vPositions[1][3]+1,vPositions[2][3]+1, vPositions[3][3]+1);
			
            this->mKMatrices.push_back(Ke);
				
            vPositions.clear();
    }
        //cout << mass << endl;
		K.shed_col(0);
		K.shed_row(0);
		

}

void Solver::constructMe(TetrahedMesh *mesh){

    arma::Mat<double> constantMatrix;
    arma::Mat<double> Me;
    arma::Mat<double> eyeMatrix = arma::eye(3,3);

    arma::Mat<double> col1 = join_rows(join_rows(join_rows(eyeMatrix*2, eyeMatrix), eyeMatrix), eyeMatrix);
    arma::Mat<double> col2 = join_rows(join_rows(join_rows(eyeMatrix, eyeMatrix*2), eyeMatrix), eyeMatrix);
    arma::Mat<double> col3 = join_rows(join_rows(join_rows(eyeMatrix, eyeMatrix), eyeMatrix*2), eyeMatrix);
    arma::Mat<double> col4 = join_rows(join_rows(join_rows(eyeMatrix, eyeMatrix), eyeMatrix), eyeMatrix*2);

    constantMatrix = join_cols(join_cols(join_cols(col1, col2), col3), col4);

    vector<arma::Mat<double> > vPositions;

    for(unsigned int k = 0; k < mesh->getNrOfTetrahedra(); k++) {

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

            arma::Mat<double> Volmat;

			Volmat = join_rows(arma::ones(4,1),join_cols(join_cols(join_cols(trans(x0),trans(x1)),trans(x2)),trans(x3)));

            double V = det(Volmat)/6.0;

            Me = density*V*0.05*constantMatrix;

            //this->tetrahedronAssemble(this->M,Me,vPositions[0][3]+1,vPositions[1][3]+1,vPositions[2][3]+1, vPositions[3][3]+1);
			//cout << vPositions[0][3] << " " << vPositions[1][3] << " " << vPositions[2][3] << " " << vPositions[3][3] << endl;
    }

    M.shed_col(0);
    M.shed_row(0);



}

void Solver::calcNewPosition(TetrahedMesh *mesh, arma::Mat<double> Fxt)
{
    //for each vertex in tetrahedron, get position
    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();

	arma::Mat<double> C = arma::eye(this->nrOfNodes*3,this->nrOfNodes*3);
	arma::Mat<double> M1 = arma::eye(this->nrOfNodes*3,this->nrOfNodes*3)*10;
	arma::Mat<double> xLocal = arma::zeros(this->nrOfNodes*3,1);

	double alpha, beta;
    alpha = 10;
    beta = 0.2;
	C = alpha*M1+beta*K;



	for (int i = 0; i < this->nrOfNodes; i++) {

		X(3*i) = mesh->mVertices->at(i).getPosition()[0];
		X(3*i+1) = mesh->mVertices->at(i).getPosition()[1];
		X(3*i+2) =  mesh->mVertices->at(i).getPosition()[2];

		xLocal(3*i) = mOriginalPos->at(i).getPosition()[0];
		xLocal(3*i+1) = mOriginalPos->at(i).getPosition()[1];
		xLocal(3*i+2) =  mOriginalPos->at(i).getPosition()[2];

	}

	if (check == false) {
		Xpre = this->X;
		check = true;
	}
	collisionForce = arma::zeros(3*this->nrOfNodes,1);

	arma::Mat<double> outerforce(3*this->nrOfNodes,1);
	outerforce = arma::zeros(3*this->nrOfNodes,1);

	arma::Mat<double> innerforce(3*this->nrOfNodes,1);
	innerforce = arma::zeros(3*this->nrOfNodes,1);


	planeCollisionDetection(X);
	arma::Mat<double> *u;
	u = new arma::Mat<double>(this->nrOfNodes*3,1);
	
	outerforce = Fxt+grav+collisionForce;

	(*u) = this->X-xLocal;
	innerforce = this->K*(*u);

    v = arma::zeros(3*this->nrOfNodes,1);
	conjugateGradient->solve(M1+C*dt+dt*dt*(this->K),&v,(M1*Vpre*1.0 - dt*(-outerforce+innerforce)));
	Vpre = v;
	


	arma::Mat<double> x = this->X+dt*v;


	//update the local positions FOR GHOST BUNNY
	arma::Mat<double> vLocal = localVpre + (grav+collisionForce)*dt;
	localVpre = vLocal;
	arma::Mat<double> newXLocal = xLocal+dt*vLocal;


	for (int i = 0; i < this->nrOfNodes; i++) {

	    arma::Mat<double> tmp, localTmp;
	    tmp << x(3*i) << x(3*i+1) << x(3*i+2);
        localTmp << newXLocal(3*i) << newXLocal(3*i+1) << newXLocal(3*i+2);

	    mesh->mVertices->at(i).setPosition(tmp);
	    mOriginalPos->at(i).setPosition(localTmp);
	}


}

void Solver::tetrahedronAssemble(arma::Mat<double> &K ,arma::Mat<double> k, int i, int j, int m, int n){

arma::Mat<double> rowPad = arma::zeros(12,1);
arma::Mat<double> colPad = arma::zeros(1,13);



k = arma::join_rows(rowPad,k);
k = arma::join_cols(colPad,k);


K(3*i-2,3*i-2) = K(3*i-2,3*i-2) + k(1,1);
K(3*i-2,3*i-1) = K(3*i-2,3*i-1) + k(1,2);
K(3*i-2,3*i) = K(3*i-2,3*i) + k(1,3);
K(3*i-2,3*j-2) = K(3*i-2,3*j-2) + k(1,4);
K(3*i-2,3*j-1) = K(3*i-2,3*j-1) + k(1,5);
K(3*i-2,3*j) = K(3*i-2,3*j) + k(1,6);
K(3*i-2,3*m-2) = K(3*i-2,3*m-2) + k(1,7);
K(3*i-2,3*m-1) = K(3*i-2,3*m-1) + k(1,8);
K(3*i-2,3*m) = K(3*i-2,3*m) + k(1,9);
K(3*i-2,3*n-2) = K(3*i-2,3*n-2) + k(1,10);
K(3*i-2,3*n-1) = K(3*i-2,3*n-1) + k(1,11);
K(3*i-2,3*n) = K(3*i-2,3*n) + k(1,12);
K(3*i-1,3*i-2) = K(3*i-1,3*i-2) + k(2,1);
K(3*i-1,3*i-1) = K(3*i-1,3*i-1) + k(2,2);
K(3*i-1,3*i) = K(3*i-1,3*i) + k(2,3);
K(3*i-1,3*j-2) = K(3*i-1,3*j-2) + k(2,4);
K(3*i-1,3*j-1) = K(3*i-1,3*j-1) + k(2,5);
K(3*i-1,3*j) = K(3*i-1,3*j) + k(2,6);
K(3*i-1,3*m-2) = K(3*i-1,3*m-2) + k(2,7);
K(3*i-1,3*m-1) = K(3*i-1,3*m-1) + k(2,8);
K(3*i-1,3*m) = K(3*i-1,3*m) + k(2,9);
K(3*i-1,3*n-2) = K(3*i-1,3*n-2) + k(2,10);
K(3*i-1,3*n-1) = K(3*i-1,3*n-1) + k(2,11);
K(3*i-1,3*n) = K(3*i-1,3*n) + k(2,12);
K(3*i,3*i-2) = K(3*i,3*i-2) + k(3,1);
K(3*i,3*i-1) = K(3*i,3*i-1) + k(3,2);
K(3*i,3*i) = K(3*i,3*i) + k(3,3);
K(3*i,3*j-2) = K(3*i,3*j-2) + k(3,4);
K(3*i,3*j-1) = K(3*i,3*j-1) + k(3,5);
K(3*i,3*j) = K(3*i,3*j) + k(3,6);
K(3*i,3*m-2) = K(3*i,3*m-2) + k(3,7);
K(3*i,3*m-1) = K(3*i,3*m-1) + k(3,8);
K(3*i,3*m) = K(3*i,3*m) + k(3,9);
K(3*i,3*n-2) = K(3*i,3*n-2) + k(3,10);
K(3*i,3*n-1) = K(3*i,3*n-1) + k(3,11);
K(3*i,3*n) = K(3*i,3*n) + k(3,12);
K(3*j-2,3*i-2) = K(3*j-2,3*i-2) + k(4,1);
K(3*j-2,3*i-1) = K(3*j-2,3*i-1) + k(4,2);
K(3*j-2,3*i) = K(3*j-2,3*i) + k(4,3);
K(3*j-2,3*j-2) = K(3*j-2,3*j-2) + k(4,4);
K(3*j-2,3*j-1) = K(3*j-2,3*j-1) + k(4,5);
K(3*j-2,3*j) = K(3*j-2,3*j) + k(4,6);
K(3*j-2,3*m-2) = K(3*j-2,3*m-2) + k(4,7);
K(3*j-2,3*m-1) = K(3*j-2,3*m-1) + k(4,8);
K(3*j-2,3*m) = K(3*j-2,3*m) + k(4,9);
K(3*j-2,3*n-2) = K(3*j-2,3*n-2) + k(4,10);
K(3*j-2,3*n-1) = K(3*j-2,3*n-1) + k(4,11);
K(3*j-2,3*n) = K(3*j-2,3*n) + k(4,12);
K(3*j-1,3*i-2) = K(3*j-1,3*i-2) + k(5,1);
K(3*j-1,3*i-1) = K(3*j-1,3*i-1) + k(5,2);
K(3*j-1,3*i) = K(3*j-1,3*i) + k(5,3);
K(3*j-1,3*j-2) = K(3*j-1,3*j-2) + k(5,4);
K(3*j-1,3*j-1) = K(3*j-1,3*j-1) + k(5,5);
K(3*j-1,3*j) = K(3*j-1,3*j) + k(5,6);
K(3*j-1,3*m-2) = K(3*j-1,3*m-2) + k(5,7);
K(3*j-1,3*m-1) = K(3*j-1,3*m-1) + k(5,8);
K(3*j-1,3*m) = K(3*j-1,3*m) + k(5,9);
K(3*j-1,3*n-2) = K(3*j-1,3*n-2) + k(5,10);
K(3*j-1,3*n-1) = K(3*j-1,3*n-1) + k(5,11);
K(3*j-1,3*n) = K(3*j-1,3*n) + k(5,12);
K(3*j,3*i-2) = K(3*j,3*i-2) + k(6,1);
K(3*j,3*i-1) = K(3*j,3*i-1) + k(6,2);
K(3*j,3*i) = K(3*j,3*i) + k(6,3);
K(3*j,3*j-2) = K(3*j,3*j-2) + k(6,4);
K(3*j,3*j-1) = K(3*j,3*j-1) + k(6,5);
K(3*j,3*j) = K(3*j,3*j) + k(6,6);
K(3*j,3*m-2) = K(3*j,3*m-2) + k(6,7);
K(3*j,3*m-1) = K(3*j,3*m-1) + k(6,8);
K(3*j,3*m) = K(3*j,3*m) + k(6,9);
K(3*j,3*n-2) = K(3*j,3*n-2) + k(6,10);
K(3*j,3*n-1) = K(3*j,3*n-1) + k(6,11);
K(3*j,3*n) = K(3*j,3*n) + k(6,12);
K(3*m-2,3*i-2) = K(3*m-2,3*i-2) + k(7,1);
K(3*m-2,3*i-1) = K(3*m-2,3*i-1) + k(7,2);
K(3*m-2,3*i) = K(3*m-2,3*i) + k(7,3);
K(3*m-2,3*j-2) = K(3*m-2,3*j-2) + k(7,4);
K(3*m-2,3*j-1) = K(3*m-2,3*j-1) + k(7,5);
K(3*m-2,3*j) = K(3*m-2,3*j) + k(7,6);
K(3*m-2,3*m-2) = K(3*m-2,3*m-2) + k(7,7);
K(3*m-2,3*m-1) = K(3*m-2,3*m-1) + k(7,8);
K(3*m-2,3*m) = K(3*m-2,3*m) + k(7,9);
K(3*m-2,3*n-2) = K(3*m-2,3*n-2) + k(7,10);
K(3*m-2,3*n-1) = K(3*m-2,3*n-1) + k(7,11);
K(3*m-2,3*n) = K(3*m-2,3*n) + k(7,12);
K(3*m-1,3*i-2) = K(3*m-1,3*i-2) + k(8,1);
K(3*m-1,3*i-1) = K(3*m-1,3*i-1) + k(8,2);
K(3*m-1,3*i) = K(3*m-1,3*i) + k(8,3);
K(3*m-1,3*j-2) = K(3*m-1,3*j-2) + k(8,4);
K(3*m-1,3*j-1) = K(3*m-1,3*j-1) + k(8,5);
K(3*m-1,3*j) = K(3*m-1,3*j) + k(8,6);
K(3*m-1,3*m-2) = K(3*m-1,3*m-2) + k(8,7);
K(3*m-1,3*m-1) = K(3*m-1,3*m-1) + k(8,8);
K(3*m-1,3*m) = K(3*m-1,3*m) + k(8,9);
K(3*m-1,3*n-2) = K(3*m-1,3*n-2) + k(8,10);
K(3*m-1,3*n-1) = K(3*m-1,3*n-1) + k(8,11);
K(3*m-1,3*n) = K(3*m-1,3*n) + k(8,12);
K(3*m,3*i-2) = K(3*m,3*i-2) + k(9,1);
K(3*m,3*i-1) = K(3*m,3*i-1) + k(9,2);
K(3*m,3*i) = K(3*m,3*i) + k(9,3);
K(3*m,3*j-2) = K(3*m,3*j-2) + k(9,4);
K(3*m,3*j-1) = K(3*m,3*j-1) + k(9,5);
K(3*m,3*j) = K(3*m,3*j) + k(9,6);
K(3*m,3*m-2) = K(3*m,3*m-2) + k(9,7);
K(3*m,3*m-1) = K(3*m,3*m-1) + k(9,8);
K(3*m,3*m) = K(3*m,3*m) + k(9,9);
K(3*m,3*n-2) = K(3*m,3*n-2) + k(9,10);
K(3*m,3*n-1) = K(3*m,3*n-1) + k(9,11);
K(3*m,3*n) = K(3*m,3*n) + k(9,12);
K(3*n-2,3*i-2) = K(3*n-2,3*i-2) + k(10,1);
K(3*n-2,3*i-1) = K(3*n-2,3*i-1) + k(10,2);
K(3*n-2,3*i) = K(3*n-2,3*i) + k(10,3);
K(3*n-2,3*j-2) = K(3*n-2,3*j-2) + k(10,4);
K(3*n-2,3*j-1) = K(3*n-2,3*j-1) + k(10,5);
K(3*n-2,3*j) = K(3*n-2,3*j) + k(10,6);
K(3*n-2,3*m-2) = K(3*n-2,3*m-2) + k(10,7);
K(3*n-2,3*m-1) = K(3*n-2,3*m-1) + k(10,8);
K(3*n-2,3*m) = K(3*n-2,3*m) + k(10,9);
K(3*n-2,3*n-2) = K(3*n-2,3*n-2) + k(10,10);
K(3*n-2,3*n-1) = K(3*n-2,3*n-1) + k(10,11);
K(3*n-2,3*n) = K(3*n-2,3*n) + k(10,12);
K(3*n-1,3*i-2) = K(3*n-1,3*i-2) + k(11,1);
K(3*n-1,3*i-1) = K(3*n-1,3*i-1) + k(11,2);
K(3*n-1,3*i) = K(3*n-1,3*i) + k(11,3);
K(3*n-1,3*j-2) = K(3*n-1,3*j-2) + k(11,4);
K(3*n-1,3*j-1) = K(3*n-1,3*j-1) + k(11,5);
K(3*n-1,3*j) = K(3*n-1,3*j) + k(11,6);
K(3*n-1,3*m-2) = K(3*n-1,3*m-2) + k(11,7);
K(3*n-1,3*m-1) = K(3*n-1,3*m-1) + k(11,8);
K(3*n-1,3*m) = K(3*n-1,3*m) + k(11,9);
K(3*n-1,3*n-2) = K(3*n-1,3*n-2) + k(11,10);
K(3*n-1,3*n-1) = K(3*n-1,3*n-1) + k(11,11);
K(3*n-1,3*n) = K(3*n-1,3*n) + k(11,12);
K(3*n,3*i-2) = K(3*n,3*i-2) + k(12,1);
K(3*n,3*i-1) = K(3*n,3*i-1) + k(12,2);
K(3*n,3*i) = K(3*n,3*i) + k(12,3);
K(3*n,3*j-2) = K(3*n,3*j-2) + k(12,4);
K(3*n,3*j-1) = K(3*n,3*j-1) + k(12,5);
K(3*n,3*j) = K(3*n,3*j) + k(12,6);
K(3*n,3*m-2) = K(3*n,3*m-2) + k(12,7);
K(3*n,3*m-1) = K(3*n,3*m-1) + k(12,8);
K(3*n,3*m) = K(3*n,3*m) + k(12,9);
K(3*n,3*n-2) = K(3*n,3*n-2) + k(12,10);
K(3*n,3*n-1) = K(3*n,3*n-1) + k(12,11);
K(3*n,3*n) = K(3*n,3*n) + k(12,12);

k.shed_col(0);
k.shed_row(0);

}

//Detects collisions
void Solver::planeCollisionDetection(arma::Mat<double> X)
{
    double planeY = 0.3;
	

    for(int i = 0; i < (X.n_rows/3.0); i++)
    {
        if(X(i*3+2) > planeY) {//Check the Y value
            planeCollisionHandler(i*3+2);
			
		}
    }
}
//Handles collision by calculating a new force
void Solver::planeCollisionHandler(unsigned int forceIndex)
{

    //F = m(v_final - v_initial) == F = -m(v_initial + v_initial) :)))

    collisionForce(forceIndex) = -15*( v(forceIndex))/dt;
	

	
    
}

arma::Mat<double> Solver::TetrahedronElementStiffness(float E,float NU,arma::Mat<double> x1,arma::Mat<double> x2, arma::Mat<double> x3, arma::Mat<double> x4) {

arma::Mat<double> xyz = arma::zeros(4,4); 
xyz << 1 << x1(0) << x1(1) << x1(2) << arma::endr 
	<< 1 << x2(0) << x2(1) << x2(2) << arma::endr 
	<< 1 << x3(0) << x3(1) << x3(2) << arma::endr 
	<< 1 << x4(0) << x4(1) << x4(2) << arma::endr; 


double V = arma::det(xyz)/6;

arma::Mat<double> mbeta1 = arma::zeros(3,3); 
mbeta1 << 1 << x2(1) << x2(2) <<  arma::endr 
	   << 1 << x3(1) << x3(2) <<  arma::endr 
	   << 1 << x4(1) << x4(2) <<  arma::endr;



arma::Mat<double> mbeta2 = arma::zeros(3,3); 
mbeta2 << 1 << x1(1) << x1(2) <<  arma::endr 
	   << 1 << x3(1) << x3(2) <<  arma::endr 
	   << 1 << x4(1) << x4(2) <<  arma::endr;

arma::Mat<double> mbeta3 = arma::zeros(3,3); 
mbeta3 << 1 << x1(1) << x1(2) <<  arma::endr 
	   << 1 << x2(1) << x2(2) <<  arma::endr 
	   << 1 << x4(1) << x4(2) <<  arma::endr;

arma::Mat<double> mbeta4 = arma::zeros(3,3); 
mbeta4 << 1 << x1(1) << x1(2) <<  arma::endr 
	   << 1 << x2(1) << x2(2) <<  arma::endr 
	   << 1 << x3(1) << x3(2) <<  arma::endr;

arma::Mat<double> mgamma1 = arma::zeros(3,3); 
mgamma1 << 1 << x2(0) << x2(2) <<  arma::endr 
	   << 1 << x3(0) << x3(2) <<  arma::endr 
	   << 1 << x4(0) << x4(2) <<  arma::endr;

arma::Mat<double> mgamma2 = arma::zeros(3,3); 
mgamma2 << 1 << x1(0) << x1(2) <<  arma::endr 
	   << 1 << x3(0) << x3(2) <<  arma::endr 
	   << 1 << x4(0) << x4(2) <<  arma::endr;

arma::Mat<double> mgamma3 = arma::zeros(3,3); 
mgamma3 << 1 << x1(0) << x1(2) <<  arma::endr 
	   << 1 << x2(0) << x2(2) <<  arma::endr 
	   << 1 << x4(0) << x4(2) <<  arma::endr;

arma::Mat<double> mgamma4 = arma::zeros(3,3); 
mgamma4 << 1 << x1(0) << x1(2) <<  arma::endr 
	   << 1 << x2(0) << x2(2) <<  arma::endr 
	   << 1 << x3(0) << x3(2) <<  arma::endr;

arma::Mat<double> mdelta1 = arma::zeros(3,3); 
mdelta1 << 1 << x2(0) << x2(1) <<  arma::endr 
	   << 1 << x3(0) << x3(1) <<  arma::endr 
	   << 1 << x4(0) << x4(1) <<  arma::endr;

arma::Mat<double> mdelta2 = arma::zeros(3,3); 
mdelta2 << 1 << x1(0) << x1(1) <<  arma::endr 
	   << 1 << x3(0) << x3(1) <<  arma::endr 
	   << 1 << x4(0) << x4(1) <<  arma::endr;

arma::Mat<double> mdelta3 = arma::zeros(3,3); 
mdelta3 << 1 << x1(0) << x1(1) <<  arma::endr 
	   << 1 << x2(0) << x2(1) <<  arma::endr 
	   << 1 << x4(0) << x4(1) <<  arma::endr;

arma::Mat<double> mdelta4 = arma::zeros(3,3); 
mdelta4 << 1 << x1(0) << x1(1) <<  arma::endr 
	   << 1 << x2(0) << x2(1) <<  arma::endr 
	   << 1 << x3(0) << x3(1) <<  arma::endr;



double beta1 = -1*det(mbeta1);

cout<< endl << mbeta4 << endl;
double beta2 = det(mbeta2);
double beta3 = -1*det(mbeta3);
double beta4 = det(mbeta4);
double gamma1 = det(mgamma1);
double gamma2 = -1*det(mgamma2);
double gamma3 = det(mgamma3);
double gamma4 = -1*det(mgamma4);
double delta1 = -1*det(mdelta1);
double delta2 = det(mdelta2);
double delta3 = -1*det(mdelta3);
double delta4 = det(mdelta4);

arma::Mat<double> B1 = arma::zeros(3,6); 
B1     << beta1 << 0 << 0 <<  arma::endr 
	   << 0 << gamma1 << 0 << arma::endr 
	   << 0 << 0 << delta1 <<  arma::endr
	   << gamma1 << beta1 << 0 << arma::endr 
	   << 0 << delta1 << gamma1 << arma::endr  
	   << delta1 << 0 << beta1 << arma::endr;
		
arma::Mat<double> B2 = arma::zeros(3,6);
B2     << beta2 << 0 << 0 <<  arma::endr 
	   << 0 << gamma2 << 0 << arma::endr 
	   << 0 << 0 << delta2 <<  arma::endr
	   << gamma2 << beta2 << 0 << arma::endr 
	   << 0 << delta2 << gamma2 << arma::endr  
	   << delta2 << 0 << beta2 << arma::endr;

arma::Mat<double> B3 = arma::zeros(3,6);
B3     << beta3 << 0 << 0 <<  arma::endr 
	   << 0 << gamma3 << 0 << arma::endr 
	   << 0 << 0 << delta3 <<  arma::endr
	   << gamma3 << beta3 << 0 << arma::endr 
	   << 0 << delta3 << gamma3 << arma::endr  
	   << delta3 << 0 << beta3 << arma::endr;

arma::Mat<double> B4 = arma::zeros(3,6);
B4     << beta4 << 0 << 0 <<  arma::endr 
	   << 0 << gamma4 << 0 << arma::endr 
	   << 0 << 0 << delta4 <<  arma::endr
	   << gamma4 << beta4 << 0 << arma::endr 
	   << 0 << delta4 << gamma4 << arma::endr  
	   << delta4 << 0 << beta4 << arma::endr;

arma::Mat<double> B = arma::zeros(6,6);

B = arma::join_rows(arma::join_rows(arma::join_rows(B1,B2),B3),B4)/(6*V);


double a = E*((1-NU)/((1+NU)*(1-2*NU)));
double b = E*(NU/((1+NU)*(1-2*NU)));
double c = E*((1-2*NU)/(2*((1+NU)*(1-2*NU))));

arma::Mat<double> AB;
arma::Mat<double> C;
AB << a << b << b << arma::endr
<< b << a << b << arma::endr
<< b << b << a << arma::endr;

C  << c << 0 << 0 << arma::endr
<< 0 << c << 0 << arma::endr
<< 0 << 0 << c << arma::endr;

arma::Mat<double> fillout = arma::zeros(3,3);
arma::Mat<double> D; 
D = arma::join_cols(arma::join_rows(AB,fillout),arma::join_rows(fillout,C));



arma::Mat<double> K =  V*trans(B)*D*B;

return K;
}
