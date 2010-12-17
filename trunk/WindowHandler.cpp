#include "WindowHandler.h"

int X;
int Y;
int g_InflictForce;
int g_NextNode;
int g_NextFourNode;
int g_AllowFracture;
int g_Gravity;


//  Callback function called when the 'AutoRotate' variable value of the tweak bar has changed
void TW_CALL SetFractureModeCB(const void *value, void *clientData)
{
    (void)clientData; // unused

    g_AllowFracture = *(const int *)(value); // copy value to g_AutoRotate

}


//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL GetFractureModeCB(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)(value) =  g_AllowFracture; // copy g_AutoRotate to value
}

//  Callback function called when the 'AutoRotate' variable value of the tweak bar has changed
void TW_CALL SetForceInflictCB(const void *value, void *clientData)
{
    (void)clientData; // unused

    g_InflictForce = *(const int *)(value); // copy value to g_AutoRotate

}


//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL GetGravityCB(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)(value) =  g_Gravity; // copy g_AutoRotate to value
}

void TW_CALL SetGravityCB(const void *value, void *clientData)
{
    (void)clientData; // unused

    g_Gravity = *(const int *)(value); // copy value to g_AutoRotate

}


//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL GetForceInflictCB(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)(value) =  g_InflictForce; // copy g_AutoRotate to value
}

//  Callback function called when the 'AutoRotate' variable value of the tweak bar has changed
void TW_CALL SetNextNodeCB(const void *value, void *clientData)
{
    (void)clientData; // unused

    g_NextNode = *(const int *)(value); // copy value to g_AutoRotate

}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL GetNextNodeCB(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)(value) =  g_NextNode; // copy g_AutoRotate to value
}

void TW_CALL SetNextFourNodeCB(const void *value, void *clientData)
{
    (void)clientData; // unused

     g_NextFourNode = *(const int *)(value); // copy value to g_AutoRotate

}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL GetNextFourNodeCB(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)(value) =  g_NextFourNode; // copy g_AutoRotate to value
}

void TW_CALL normalCB(void *clientdata)
{
	WindowHandler* test = static_cast<WindowHandler*>(clientdata);
	test->volumeGenerator->changeNormalRenderMode();
}

void TW_CALL edgeCB(void *clientdata)
{
	WindowHandler* test = static_cast<WindowHandler*>(clientdata);
	test->volumeGenerator->changeEdgeRenderMode();
}

void TW_CALL triangleCB(void *clientdata)
{
	WindowHandler* test = static_cast<WindowHandler*>(clientdata);
	test->volumeGenerator->changeTriangleRenderMode();
}





WindowHandler::WindowHandler(void)
{

    windowWidth = 800;
    windowHeight = 700;
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("");

	bar = TwNewBar("TweakBar");
	TwWindowSize(800, 700);
	TwDefine(" TweakBar size='230 600' color='0 0 0' "); // change default tweak bar size and color

	this->g_Rotation[0] =  0.0f;
	this->g_Rotation[1] =  0.0f;
	this->g_Rotation[2] =  0.0f;
	this->g_Rotation[3] =  1.0f;

	g_InflictForce = 0;
	g_NormalMode = 0;
	g_EdgeMode = 0;
	g_TriangleMode = 0;

	g_NextNode = 0;
	g_Dampening = 0;
	g_Mass = 23;
	g_E = 2300000;
	g_vn = 0.3;


	this->g_RotateStart[0] = 0.0f;
	this->g_RotateStart[1] = 0.0f;
	this->g_RotateStart[2] = 0.0f;
	this->g_RotateStart[3] = 1.0f;

	this->g_ForceDirection[0] = 0.5f;
	this->g_ForceDirection[1] = 0.5f;
	this->g_ForceDirection[2] = 0.5f;
	this->g_Force = 0.0f;
	this->g_fractureThresh = 20000;
	this->g_alpha = 10;
	this->g_beta = 0.06;


    //Set arcball
    eye.setVec( 0.0f, 0.1f, 0.5f );
    center.setVec( 0.0f, 0.0f, 0.0f );
    up.setVec( 0.0f, 1.0f, 0.0f );

    SPHERE_RADIUS = 1.0f;
    PI = 3.141592654f;
    buttonPressed = -1;



	// Add callback to toggle auto-rotate mode (callback functions are defined above).
    TwAddVarCB(bar, "Inflict Force", TW_TYPE_BOOL32, SetForceInflictCB, GetForceInflictCB, NULL,
               " label='Inflict Force' key=space help='Toggle Force mode.' ");

	   // Add 'g_Zoom' to 'bar': this is a modifable (RW) variable of type TW_TYPE_FLOAT. Its key shortcuts are [z] and [Z].
    TwAddVarRW(bar, "Force", TW_TYPE_FLOAT, &this->g_Force,
               " min=0.00 max=4000 ; step=1.0 keyIncr=z keyDecr=Z help='Force applied on Node' ");

	TwAddVarRW(bar, "Force Direction", TW_TYPE_DIR3F, &this->g_ForceDirection,
               " label='Force direction' open help='Change Force Direction' ");

               	// Add callback to toggle auto-rotate mode (callback functions are defined above).
    TwAddVarCB(bar, "Next Node", TW_TYPE_BOOL32, SetNextNodeCB, GetNextNodeCB, NULL,
               " label='Next Node' key=space help='Toggle Next Node.' ");


    TwAddVarCB(bar, "Gravity", TW_TYPE_BOOL32, SetGravityCB, GetGravityCB, NULL,
               " label='Gravity' key=space help='Toggle Force mode.' ");

               // Add callback to toggle auto-rotate mode (callback functions are defined above).
    TwAddVarCB(bar, "Next Four Nodes", TW_TYPE_BOOL32, SetNextFourNodeCB, GetNextFourNodeCB, NULL,
               " label='Attach' key=space help='Toggle Next Node.' ");

    // Add 'g_Rotation' to 'bar': this is a variable of type TW_TYPE_QUAT4F which defines the object's orientation
    TwAddVarRW(bar, "ObjRotation", TW_TYPE_QUAT4F, &this->g_Rotation,
               " label='Object rotation' open help='Change the object orientation.' ");

	TwAddSeparator(bar, NULL, " group='MaterialSettings' ");
	TwAddVarRW(bar, "Mass", TW_TYPE_FLOAT, &this->g_Mass,
               "group='MaterialSettings' min=0.1 max=4000;");
	TwAddVarRW(bar, "Alpha", TW_TYPE_FLOAT, &this->g_alpha,
               "group='MaterialSettings' min=-100.00 max=100 step=0.1;");
	TwAddVarRW(bar, "Beta ", TW_TYPE_FLOAT, &this->g_beta,
               "group='MaterialSettings' min=0.00 max=100 step = 0.001;");
	TwAddVarRW(bar, "Young's modulus ", TW_TYPE_FLOAT, &this->g_E,
               "group='MaterialSettings' min=1000 max=30000000000 step = 1000;");
	TwAddVarRW(bar, "Possion ratio ", TW_TYPE_FLOAT, &this->g_vn,
               "group='MaterialSettings' min=0.0 max=0.49 step = 0.01;");

	TwAddSeparator(bar, NULL, " group='RenderSettings' ");
	TwAddButton(bar, "Normal Mode", normalCB, this, " group='RenderSettings' ");
	TwAddButton(bar, "Edge Mode", edgeCB, this, " group='RenderSettings' ");
	TwAddButton(bar, "Triangle Mode", triangleCB, this, " group='RenderSettings' ");

	TwAddSeparator(bar, NULL, " group='FractureSettings' ");
		// Add callback to toggle auto-rotate mode (callback functions are defined above).
    TwAddVarCB(bar, "Fracture Mode", TW_TYPE_BOOL32, SetFractureModeCB, GetFractureModeCB, NULL,
               " label='Inflict Force' key=space help='Toggle Force mode.' ");

	TwAddVarRW(bar, "Fracture limit", TW_TYPE_FLOAT, &this->g_fractureThresh,
               "group='FractureSettings' min=10000 max=40000;");


}
WindowHandler::~WindowHandler(void)
{

}

void WindowHandler::display()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    //RenderFirstPass(); //Deformation Simulation
    RenderSecondPass(); //Render the actual graphics
	TwDraw();
    glutSwapBuffers();

}
//Function for the deformation simulation
void WindowHandler::RenderFirstPass()
{
   	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
        glLoadIdentity();
        glViewport(0, 0, textureSize, textureSize);
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
            glLoadIdentity();
            gluOrtho2D(0.0, textureSize, 0.0, textureSize);
            glEnable(GL_TEXTURE_2D);

            glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
            GLenum drawBuffers[] = {GL_COLOR_ATTACHMENT0_EXT,GL_COLOR_ATTACHMENT1_EXT};
            glDrawBuffers(2, drawBuffers);

            deformationShader->use();
            deformationShader->sendUniformTexture("positionTex",0);

            //Quad for the shader to use
            drawQuad();
            deformationShader->disable();

            glDisable(GL_TEXTURE_2D);
            glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glViewport(0, 0, windowWidth,windowHeight);
	glFlush();

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

}
//Function for the rendering to the screen
void WindowHandler::RenderSecondPass()
{

  float mat[4*4];
  glEnable(GL_CULL_FACE);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);

  keyHandler();
  this->solver->setParameter(this->g_Mass, this->g_fractureThresh, g_AllowFracture, g_alpha, g_beta, g_E, g_vn);

  glPushMatrix();

  //arma::Mat<double> temp3;
//  temp3 = temp*temp2;


  ConvertQuaternionToMatrix(g_Rotation, mat);


  glMultMatrixf(mat);

  solver->calcNewPosition(volumeGenerator->getTetrahedMesh(), this->Fxt);
  this->Fxt = arma::zeros(this->Fxt.n_rows,this->Fxt.n_cols);

  lightShader->use();
  volumeGenerator->render();
  lightShader->disable();

  //

  this->drawForceArrow();
  glPopMatrix();

}

void WindowHandler::setupTextures()
{
    ///Example, if we need the position of the vertices
    //Get the position data for each vertex
	/*GLfloat *positionData = volumeGenerator->getTetrahedMesh().GetVertexArray();


    //Create the position texture that will be sent to the shader for integration
	glGenTextures(1, &positionTexID);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, positionTexID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, textureSize, textureSize, 0, GL_RGBA, GL_FLOAT, &positionData[0]);

    //Done with the position data
    delete [] positionData;

    //generate a framebuffer object and bind the textures to it.
	glGenFramebuffersEXT(1, &fbo);
  	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	glBindTexture(GL_TEXTURE_2D, positionTexID);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, positionTexID, 0);
    ///END EXAMPLE

	if (glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT) != GL_FRAMEBUFFER_COMPLETE_EXT)
	printf("ERROR - Incomplete FrameBuffer\n");

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);*/
}

void WindowHandler::init()
{

    lightShader = new Shader();
    lightShader->load("Shader/vertexPhongShader.glsl","Shader/fragmentPhongShader.glsl");

    deformationShader = new Shader();
    deformationShader->load("Shader/vertexDeformationShader.glsl","Shader/fragmentDeformationShader.glsl");

	volumeGenerator = new VolumeGenerator();
	volumeGenerator->generateVolume();


	solver = new Solver(volumeGenerator->getTetrahedMesh()->getNrOfNodes());
	Fxt = arma::zeros(volumeGenerator->getTetrahedMesh()->getNrOfNodes()*3,1);

	//volumeGenerator->subdivide();

	//For the deformation
    textureSize = volumeGenerator->getTetrahedMesh()->GetVertexArraySize(); //a texture is optimal if 2^n large
    //
	this->nrOfVertices = textureSize*textureSize;
    //Setup textures used for the deformation shader
	//setupTextures();
	solver->constructKe(volumeGenerator->getTetrahedMesh());
	solver->constructMe(volumeGenerator->getTetrahedMesh());


}

// Routine to set a quaternion from a rotation axis and angle
// ( input axis = float[3] angle = float  output: quat = float[4] )
void WindowHandler::SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat)
{
    float sina2, norm;
    sina2 = (float)sin(0.5f * angle);
    norm = (float)sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    quat[0] = sina2 * axis[0] / norm;
    quat[1] = sina2 * axis[1] / norm;
    quat[2] = sina2 * axis[2] / norm;
    quat[3] = (float)cos(0.5f * angle);

}


// Routine to convert a quaternion to a 4x4 matrix
// ( input: quat = float[4]  output: mat = float[4*4] )
void WindowHandler::ConvertQuaternionToMatrix(const float *quat, float *mat)
{
    float yy2 = 2.0f * quat[1] * quat[1];
    float xy2 = 2.0f * quat[0] * quat[1];
    float xz2 = 2.0f * quat[0] * quat[2];
    float yz2 = 2.0f * quat[1] * quat[2];
    float zz2 = 2.0f * quat[2] * quat[2];
    float wz2 = 2.0f * quat[3] * quat[2];
    float wy2 = 2.0f * quat[3] * quat[1];
    float wx2 = 2.0f * quat[3] * quat[0];
    float xx2 = 2.0f * quat[0] * quat[0];
    mat[0*4+0] = - yy2 - zz2 + 1.0f;
    mat[0*4+1] = xy2 + wz2;
    mat[0*4+2] = xz2 - wy2;
    mat[0*4+3] = 0;
    mat[1*4+0] = xy2 - wz2;
    mat[1*4+1] = - xx2 - zz2 + 1.0f;
    mat[1*4+2] = yz2 + wx2;
    mat[1*4+3] = 0;
    mat[2*4+0] = xz2 + wy2;
    mat[2*4+1] = yz2 - wx2;
    mat[2*4+2] = - xx2 - yy2 + 1.0f;
    mat[2*4+3] = 0;
    mat[3*4+0] = mat[3*4+1] = mat[3*4+2] = 0;
    mat[3*4+3] = 1;
}


// Routine to multiply 2 quaternions (ie, compose rotations)
// ( input q1 = float[4] q2 = float[4]  output: qout = float[4] )
void WindowHandler::MultiplyQuaternions(const float *q1, const float *q2, float *qout)
{
    float qr[4];
	qr[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
	qr[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
	qr[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
	qr[3]  = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
    qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}

void WindowHandler::reshape(int w, int h)
{
    windowHeight=h;
    windowWidth=w;
    glClearColor(0.0, 0.0, 0.0, 1.0);

    float aspect_ratio = (float) windowWidth / (float) windowHeight;

    glViewport(0, 0, windowWidth, windowHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( 60.0f, aspect_ratio, 0.05f, 650.0f );

    arcball_setzoom( SPHERE_RADIUS, eye, up );

    glMatrixMode(GL_MODELVIEW);
    //glLoadIdentity() ;
	gluLookAt(
        eye.x, eye.y, eye.z,   //eye
        center.x, center.y, center.z,   //lookat
        up.x, up.y, up.z );  //up vector

}
void WindowHandler::idle()
{
    showFPS();
    glutPostRedisplay();
}

void WindowHandler::mouseButtonEvent(int button, int state, int x, int y)
{
    buttonPressed = button;

	if( state == GLUT_DOWN ) {

		if( button == GLUT_LEFT_BUTTON) {
			int invert_y = (windowHeight - y) - 1; // OpenGL viewport coordinates are Cartesian
			arcball_start(x,invert_y);
		}
	}

	glutPostRedisplay();
}
void WindowHandler::mouseMoveEvent(int x, int y)
{
	X = x;
	Y = y;
    if( buttonPressed == GLUT_LEFT_BUTTON || buttonPressed == GLUT_RIGHT_BUTTON || buttonPressed == GLUT_MIDDLE_BUTTON) {

		if( buttonPressed == GLUT_LEFT_BUTTON) {
			int invert_y = (windowHeight - y) - 1;
			arcball_move(x,invert_y);
		}

		glutPostRedisplay();

	}
}

void WindowHandler::processNormalKeys(unsigned char,int,int) {


}

void WindowHandler::keyHandler() {

	/*
	if(key == 110) {

		volumeGenerator->changeNormalRenderMode();
	}

	if(key == 101) {

		volumeGenerator->changeEdgeRenderMode();
	}

	if(key == 116) {

		volumeGenerator->changeTriangleRenderMode();
	}

	if(key == 116) {

		volumeGenerator->changeTriangleRenderMode();
	}
	*/
	if (g_NextFourNode == 1) {

	this->volumeGenerator->getTetrahedMesh()->pickNextFourNode();
    g_NextFourNode = 0;
	}

	if (g_NextFourNode == 0) {

	//this->volumeGenerator->getTetrahedMesh()->pickedNodes.clear();

	}

	if (g_NextNode == 1) {

	this->volumeGenerator->getTetrahedMesh()->pickNextNode();
	g_NextNode = 0;
	}

	if (g_InflictForce == 1)
	{
		unsigned int cNode = this->volumeGenerator->getTetrahedMesh()->getCurrentNode();
		this->Fxt(cNode*3) = g_ForceDirection[0]*this->g_Force;
		this->Fxt(cNode*3+1) = g_ForceDirection[1]*this->g_Force;
		this->Fxt(cNode*3+2) = g_ForceDirection[2]*this->g_Force;
	}

		if (g_Gravity == 1)
	{
		vector<unsigned int> templist = this->volumeGenerator->getTetrahedMesh()->pickedNodes;
		unsigned int n1,n2,n3,n4;
		if (templist.size() > 3) {

		n1 = templist.at(0);
        n2 = templist.at(1);
        n3 = templist.at(2);
        n4 = templist.at(3);
        } else {

        n1 = -1;
        n2 = -1;
        n3 = -1;
        n4 = -1;
        }

            for(int j = 0; j < volumeGenerator->getTetrahedMesh()->getNrOfNodes(); j++ ) {

                if (j != n1 && j != n2 && j != n3 && j != n4){
                this->Fxt(j*3+1) = -9.82;
                }
            }

        }


	/*
	if (key == 112)
	{
	this->volumeGenerator->getTetrahedMesh()->pickNextNode();
	}
	*/


}

//Quad for the deformation texture
void WindowHandler::drawQuad()
{
    glBegin(GL_QUADS);
        glTexCoord2f(0.0, 0.0);
        glVertex2f(0.0, 0.0);
        glTexCoord2f(1.0, 0.0);
        glVertex2f( textureSize, 0.0);
        glTexCoord2f(1.0, 1.0);
        glVertex2f( textureSize,  textureSize);
        glTexCoord2f(0.0, 1.0);
        glVertex2f(0.0, textureSize);
    glEnd();
}

void WindowHandler::drawForceArrow() {


	arma::Mat<double> curNode = this->volumeGenerator->getTetrahedMesh()->pickNode();
    glColor3f(0.0,0.0,0.6);
	glPushMatrix();
	glTranslatef(curNode(0), curNode(1), curNode(2));
	glutSolidSphere(0.002,40,40);
	glPopMatrix();

	vector<unsigned int> templist = this->volumeGenerator->getTetrahedMesh()->pickedNodes;

    glColor3f(0.6,0.0,0.0);

    for(int i = 0; i < templist.size();i++){

        unsigned int currentVertex = templist.at(i);

        arma::Mat<double> temp = this->volumeGenerator->getTetrahedMesh()->mVertices->at(currentVertex).getPosition();
        glPushMatrix();
        glTranslatef(temp(0), temp(1), temp(2));
        glutSolidSphere(0.002,40,40);
        glPopMatrix();
    }

}

void WindowHandler::showFPS() {

    float t;
    float fps;

    // Get current time
    t =  (float) (glutGet(GLUT_ELAPSED_TIME)/1000.0f); // Gets number of seconds since glfwInit()
    // If one second has passed, or if this is the very first frame
    if( (t-t0) > 1.0 || frames == 0 )
    {
        fps = (float)(frames / (t-t0));
        sprintf(titlestring, "Deformable Bodies (%.1f FPS)", fps);
        glutSetWindowTitle(titlestring);
        t0 = t;
        frames = 0;
    }
    frames ++;
}
