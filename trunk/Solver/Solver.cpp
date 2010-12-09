#include "Solver.h"


Solver::Solver(int nrOfNodes){

	this->nrOfNodes = nrOfNodes;
	conjugateGradient = new ConjugateGradient(100,0.001);
	K = arma::zeros(this->nrOfNodes*3+1,this->nrOfNodes*3+1);
	M = arma::zeros(this->nrOfNodes*3+1,this->nrOfNodes*3+1);
	Xpre = arma::zeros(this->nrOfNodes*3,1);
	Vpre = arma::zeros(this->nrOfNodes*3,1);
	localVpre = arma::zeros(this->nrOfNodes*3,1);
	stress = arma::zeros(12,1);
	allowFracture = false;
	
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

	this->vn = 0.33;
    this->E = 23000.300000;




}

void Solver::setParameter(float mass, float fracture) {

	this->mass = mass;
	this->alpha = 10;
    this->beta = 0.2;
	this->FractureThresh = fracture;

}

Solver::~Solver(){

	
}


void Solver::constructKe(TetrahedMesh *mesh){

    mOriginalPos = new vector<Vertex>((*mesh->mVertices));
	arma::Mat<double> R = findRotation(mesh,0);
	
    //for each vertex in tetraheder, get positon
    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();

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

			arma::Mat<double> Ke(12,12);
			Ke = arma::zeros<arma::mat>(12,12);

			
			Ke = TetrahedronElementStiffness(E,vn, x0, x1, x2, x3);
			arma::Mat<double> X = join_cols(join_cols(join_cols(x0, x1), x2), x3);
		
			this->tetrahedronAssemble(this->K,Ke,vPositions[0][3]+1,vPositions[1][3]+1,vPositions[2][3]+1, vPositions[3][3]+1);
			
            this->mKMatrices.push_back(Ke);
				
            vPositions.clear();
    }
       
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
	
	arma::Mat<double> xLocal = arma::zeros(this->nrOfNodes*3,1);


	//Get all updated positions,
	for (int i = 0; i < this->nrOfNodes; i++) {

		X(3*i) = mesh->mVertices->at(i).getPosition()[0];
		X(3*i+1) = mesh->mVertices->at(i).getPosition()[1];
		X(3*i+2) =  mesh->mVertices->at(i).getPosition()[2];

		xLocal(3*i) = mOriginalPos->at(i).getPosition()[0];
		xLocal(3*i+1) = mOriginalPos->at(i).getPosition()[1];
		xLocal(3*i+2) =  mOriginalPos->at(i).getPosition()[2];

	}

    unsigned int nrOfTetraheds = mesh->getNrOfTetrahedra();

	vector<arma::Mat<double> > vPositions;
	this->K = arma::zeros(this->nrOfNodes*3+1,this->nrOfNodes*3+1);
	arma::Mat<double> Kf0 = arma::zeros(this->nrOfNodes*3+1,this->nrOfNodes*3+1); 
	for(int i = 0; i<nrOfTetraheds;i++){

		vPositions = mesh->getVertexPosition(i);
		arma::Mat<double> Ke = this->mKMatrices.at(i);
		arma::Mat<double> Re = arma::zeros<arma::mat>(12,12);
		Re = findRotation(mesh,i);
		
		this->tetrahedronAssemble(this->K,Re*Ke*trans(Re),vPositions[0][3]+1,vPositions[1][3]+1,vPositions[2][3]+1, vPositions[3][3]+1);
		this->tetrahedronAssemble(Kf0,Re*Ke,vPositions[0][3]+1,vPositions[1][3]+1,vPositions[2][3]+1, vPositions[3][3]+1);
		if (allowFracture == true) {
		arma::Mat<double> stresstensor = calculateStress(mesh,i, this->E,this->vn, Re);
		crackIT(mesh, i, stresstensor);
		}
	
	}
		Kf0.shed_col(0);
		Kf0.shed_row(0);
		this->K.shed_col(0);
		this->K.shed_row(0);

	//Init mass and dampening matrix
	this->M = arma::eye(this->nrOfNodes*3,this->nrOfNodes*3)*this->mass;
	this->C = arma::eye(this->nrOfNodes*3,this->nrOfNodes*3);
	C = this->alpha*M+this->beta*K;



	collisionForce = arma::zeros(3*this->nrOfNodes,1);

	arma::Mat<double> outerforce(3*this->nrOfNodes,1);
	outerforce = arma::zeros(3*this->nrOfNodes,1);

	arma::Mat<double> innerforce(3*this->nrOfNodes,1);
	innerforce = arma::zeros(3*this->nrOfNodes,1);


	planeCollisionDetection(X);
	arma::Mat<double> *u;
	u = new arma::Mat<double>(this->nrOfNodes*3,1);
	
	outerforce = Fxt+grav+collisionForce;

	
	innerforce =  this->K*this->X-Kf0*xLocal;

    v = arma::zeros(3*this->nrOfNodes,1);
	conjugateGradient->solve(this->M+this->C*dt+dt*dt*(this->K),&v,(this->M*Vpre*1.0 - dt*(-outerforce+innerforce)));
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

void Solver::changeFracture() {

	this->allowFracture = true;
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

arma::Mat<double> B = arma::zeros(6,6);

B = calculateB(x1, x2, x3,  x4);

mBMatrices.push_back(B);

arma::Mat<double> D; 

D = calculateD(E,NU);

arma::Mat<double> K =  V*trans(B)*D*B;

return K;
}

arma::Mat<double> Solver::calculateB(arma::Mat<double> x1,arma::Mat<double> x2, arma::Mat<double> x3, arma::Mat<double> x4) {

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

return B;

}

arma::Mat<double> Solver::calculateD(float E,float NU) {

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

	return D;
}

arma::Mat<double> Solver::findRotation(TetrahedMesh *mesh,unsigned int theInd) {

	vector<arma::Mat<double>> vDisPos;
	vector<arma::Mat<double>> vPos;
	arma::Mat<double> X, P;
	vDisPos = mesh->getVertexPosition(0);
	

	
            arma::Col<double> p0(3);
            arma::Col<double> p1(3);
            arma::Col<double> p2(3);
			arma::Col<double> p3(3);

			arma::Col<double> x0(3);
            arma::Col<double> x1(3);
            arma::Col<double> x2(3);
			arma::Col<double> x3(3);

            p0(0) = vDisPos[0][0];
            p0(1) = vDisPos[0][1];
            p0(2) = vDisPos[0][2];

            p1(0) = vDisPos[1][0];
            p1(1) = vDisPos[1][1];
            p1(2) = vDisPos[1][2];

            p2(0) = vDisPos[2][0];
            p2(1) = vDisPos[2][1];
            p2(2) = vDisPos[2][2];

            p3(0) = vDisPos[3][0];
            p3(1) = vDisPos[3][1];
            p3(2) = vDisPos[3][2];

			
			x0 = trans(this->mOriginalPos->at(vDisPos[0][3]).getPosition());
			x1 = trans(this->mOriginalPos->at(vDisPos[1][3]).getPosition());
			x2 = trans(this->mOriginalPos->at(vDisPos[2][3]).getPosition());
			x3 = trans(this->mOriginalPos->at(vDisPos[3][3]).getPosition());
			arma::Row<double> pad = arma::ones(1,4);
			P = join_rows(join_rows((p1-p0),(p2-p0)), (p3-p0));
			X = join_rows(join_rows((x1-x0),(x2-x0)), (x3-x0));
			
			arma::Mat<double> R = calculateRotation(X, P);
			arma::Mat<double> Z = arma::zeros(3,3);

			arma::Mat<double> col1 = join_rows(join_rows(join_rows(R, Z), Z), Z);
			arma::Mat<double> col2 = join_rows(join_rows(join_rows(Z, R), Z), Z);
			arma::Mat<double> col3 = join_rows(join_rows(join_rows(Z, Z), R), Z);
			arma::Mat<double> col4 = join_rows(join_rows(join_rows(Z, Z), Z), R);

			arma::Mat<double> Re = join_cols(join_cols(join_cols(col1, col2), col3), col4);
			
			//cout << endl << Re << endl;
			return Re;
			
}


arma::Mat<double> Solver::calculateRotation(arma::Mat<double> X, arma::Mat<double> P) {

	arma::Mat<double> R, A;
	A = P * inv(X);
	
	arma::Mat<double> a0 = A.col(0);
	arma::Mat<double> a1= A.col(1);
	arma::Mat<double> a2 = A.col(2);

	arma::Mat<double> r0 = a0/arma::norm(a0,2);
	arma::Mat<double> r1 = (a1-arma::dot(r0,a1))/arma::norm(a1-arma::dot(r0,a1),2);
	arma::Mat<double> r2 = arma::cross(r0,r1);

	R = join_rows(join_rows(r0,r1),r2);


return R;
}

arma::Mat<double> Solver::calculateStress(TetrahedMesh *mesh,unsigned int theInd, float E,float NU, arma::Mat<double> R) {

	vector<arma::Mat<double>> vDisPos;
	vector<arma::Mat<double>> vPos;
	arma::Mat<double> X, P;
	vDisPos = mesh->getVertexPosition(0);
	
            arma::Col<double> p0(3);
            arma::Col<double> p1(3);
            arma::Col<double> p2(3);
			arma::Col<double> p3(3);

			arma::Col<double> x0(3);
            arma::Col<double> x1(3);
            arma::Col<double> x2(3);
			arma::Col<double> x3(3);

            p0(0) = vDisPos[0][0];
            p0(1) = vDisPos[0][1];
            p0(2) = vDisPos[0][2];

            p1(0) = vDisPos[1][0];
            p1(1) = vDisPos[1][1];
            p1(2) = vDisPos[1][2];

            p2(0) = vDisPos[2][0];
            p2(1) = vDisPos[2][1];
            p2(2) = vDisPos[2][2];

            p3(0) = vDisPos[3][0];
            p3(1) = vDisPos[3][1];
            p3(2) = vDisPos[3][2];

			x0 = trans(this->mOriginalPos->at(vDisPos[0][3]).getPosition());
			x1 = trans(this->mOriginalPos->at(vDisPos[1][3]).getPosition());
			x2 = trans(this->mOriginalPos->at(vDisPos[2][3]).getPosition());
			x3 = trans(this->mOriginalPos->at(vDisPos[3][3]).getPosition());

			X = join_cols(join_cols(join_cols(x0, x1), x2), x3);
			P = join_cols(join_cols(join_cols(p0, p1), p2), p3);

			arma::Mat<double> B = this->mBMatrices.at(theInd);
			arma::Mat<double> D = this->calculateD(E,NU);

			arma::Mat<double> stress = D*B*(trans(R)*P-X);

			arma::Mat<double> stresstensor(3,3);
			stresstensor(0,0) = stress(0);
			stresstensor(1,1) = stress(1);
			stresstensor(2,2) = stress(2);

			stresstensor(0,1) = stress(3);
			stresstensor(1,0) = stress(3);

			stresstensor(1,2) = stress(4);
			stresstensor(2,1) = stress(4);

			stresstensor(0,2) = stress(5);
			stresstensor(2,0) = stress(5);



			return stresstensor;

}

arma::Mat<double> Solver::calculateLargestEIG(arma::Mat<double> stresstensor){

	arma::Col<double> eigval;
	arma::Mat<double> eigvec;
	arma::eig_sym(eigval,eigvec,stresstensor);
	
	double max = arma::max<double>(eigval);
	
	if (max > 10000.0)
	{
		arma::Col<UINT32> q2 = find(eigval ==  max,1,"first");
		
		
		return eigvec.col(q2(0));
	} 
	else {
	
		arma::Mat<double> null(1,1);
		null(0) =-1000.0;

	return null;
	}
	
}

arma::Mat<double> Solver::crackIT(TetrahedMesh *mesh, unsigned int theInd, arma::Mat<double> stresstensor) {

		arma::Mat<double> maxEig = calculateLargestEIG(stresstensor);

		if (maxEig(0) != -1000.0){

			mesh->crackStructure(theInd,maxEig);
		
		}

		arma::Mat<double> butcrap;
		return butcrap;
}