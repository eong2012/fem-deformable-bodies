#include "Solver.h"


Solver::Solver(int nrOfNodes){

	this->nrOfNodes = nrOfNodes;
	conjugateGradient = new ConjugateGradient(20.0,0.0001);
	K = arma::zeros(this->nrOfNodes*3+1,this->nrOfNodes*3+1);
	Xpre = arma::zeros(this->nrOfNodes*3,1);
	Vpre = arma::zeros(this->nrOfNodes*3,1);
	localVpre = arma::zeros(this->nrOfNodes*3,1);
	check = false;
	update = false;
	grav = arma::zeros(this->nrOfNodes*3,1);
	int i = 0;
	X = arma::zeros(this->nrOfNodes*3,1);
	dt = 0.01;
	mass = 1.1;
	this->v = arma::zeros(this->nrOfNodes*3,1);

	collisionForce = arma::zeros(this->nrOfNodes*3,1);

	while(i*3+1 < grav.n_rows) {

	grav(i*3+1) = -9.82;
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

void Solver::contstructKe(TetrahedMesh *mesh){

    mOriginalPos = new vector<Vertex>((*mesh->mVertices));

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

            cout << "--------" <<"TETRA" << k << "--------" << endl;
            cout << "vertex 1" << endl;
            cout << " x: " <<x0(0) << " y: " <<x0(1) << " z: " <<x0(2)<< endl;

            cout << "vertex 2" << endl;
            cout << " x: " <<x1(0) << " y: " <<x1(1) << " z: " <<x1(2)<< endl;

            cout << "vertex 3" << endl;
            cout << " x: " <<x2(0) << " y: " <<x2(1) << " z: " <<x2(2)<< endl;

            cout << "vertex 4" << endl;
            cout << " x: " <<x3(0) << " y: " <<x3(1) << " z: " <<x3(2)<< endl;

			arma::Mat<double> xOrgin = join_cols(join_cols(join_cols(x0,x1), x2),x3);
			this->xOrgin.push_back(xOrgin);

            arma::Mat<double> X, Volmat;

            X = join_rows(join_rows((x1-x0),(x2-x0)), (x3-x0));


			Volmat = join_rows(arma::ones(4,1),join_cols(join_cols(join_cols(trans(x0),trans(x1)),trans(x2)),trans(x3)));

            double V = det(Volmat)/6.0;
            cout << "VOLUME: " << V << endl;

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
            vn = 0.08;
            E = 0.5;

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
			//cout << vPositions[0][3] << " " << vPositions[1][3] << " " << vPositions[2][3] << " " << vPositions[3][3] << endl;
			this->tetrahedronAssemble(this->K,Ke,vPositions[0][3]+1,vPositions[1][3]+1,vPositions[2][3]+1, vPositions[3][3]+1);
            this->mKMatrices.push_back(Ke);

            vPositions.clear();
    }


		cout << "------" << endl;

		K.shed_col(0);
		K.shed_row(0);
		//cout << this->nrOfNodes  << endl;
	//cout << this->K << endl;

}

void Solver::calcNewPosition(TetrahedMesh *mesh, arma::Mat<double> Fxt)
{
     //for each vertex in tetraheder, get positon
    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();

	arma::Mat<double> M = arma::eye(this->nrOfNodes*3,this->nrOfNodes*3);
	arma::Mat<double> C = arma::eye(this->nrOfNodes*3,this->nrOfNodes*3)*1;
	arma::Mat<double> xLocal = arma::zeros(this->nrOfNodes*3,1);



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

	//NY
	arma::Mat<double> outerforce(3*this->nrOfNodes,1);
	outerforce = arma::zeros(3*this->nrOfNodes,1);

	arma::Mat<double> innerforce(3*this->nrOfNodes,1);
	innerforce = arma::zeros(3*this->nrOfNodes,1);


	planeCollisionDetection(X);
	arma::Mat<double> *u;
	u = new arma::Mat<double>(this->nrOfNodes*3,1);

	outerforce = Fxt+grav+collisionForce;

	//conjugateGradient->solve(K,u,outerforce);
	(*u) = this->X-xLocal;
	innerforce = this->K*(*u);

	cout << (*u) << endl;
	//NY

	//NY
	conjugateGradient->solve(M+C*dt+dt*dt*(this->K),&v,(M*Vpre + dt*(outerforce-innerforce)));
	Vpre = v;
	//cout << v << endl;

	//GAMLA MED EN KRAFT!!!
	//this->v = Vpre+dt*M*(Fxt+this->K*(*u)-C*Vpre);

	//cout << v << endl;


	arma::Mat<double> x = this->X+dt*v;
	//Xpre = this->X;

	//update the local positions
	arma::Mat<double> vLocal = localVpre + (grav+collisionForce)*dt;
	localVpre = vLocal;
	arma::Mat<double> newXLocal = xLocal+dt*vLocal;


	for (int i = 0; i < this->nrOfNodes; i++) {

	    arma::Mat<double> tmp, localTmp;
	    tmp << x(3*i) << x(3*i+1) << x(3*i+2);
        localTmp << newXLocal(3*i) << newXLocal(3*i+1) << newXLocal(3*i+2);

	    mesh->mVertices->at(i).setPosition(tmp);
	    mOriginalPos->at(i).setPosition(localTmp);

		//mesh->mVertices->at(i).setPosition(arma::Mat<double>(x(3*i),x(3*i+1), x(3*i+2)));

	}


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

//Detects collisions
void Solver::planeCollisionDetection(arma::Mat<double> X)
{
    double planeY = -0.7;

    for(int i = 0; i < (X.n_rows/3.0); i++)
    {
        if(X(i*3+1) < planeY) //Check the Y value
            planeCollisionHandler(i*3+1);

    }
}
//Handles collision by calculating a new force
void Solver::planeCollisionHandler(unsigned int forceIndex)
{

    //F = m(v_final - v_initial) == F = -m(v_initial + v_initial) :)))
    collisionForce(forceIndex) =-mass*(v(forceIndex) + v(forceIndex))/dt;
	Xpre = this->X;
}

