#include "WindowHandler.h"

WindowHandler::WindowHandler(void)
{

    windowWidth = 600;
    windowHeight = 600;
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("");

    //Set arcball
    eye.setVec( 0.0f, 0.0f, 3.0f );
    center.setVec( 0.0f, 0.0f, 0.0f );
    up.setVec( 0.0f, 1.0f, 0.0f );

    SPHERE_RADIUS = 1.0f;
    PI = 3.141592654f;
    buttonPressed = -1;
}
WindowHandler::~WindowHandler(void)
{

}

void WindowHandler::display()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    RenderFirstPass(); //Deformation Simulation
    RenderSecondPass(); //Render the actual graphics
    //cout << "in display" << endl;

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

  glEnable(GL_CULL_FACE);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  ///New --- Read the from the texture and store in an array
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, positionTexID);

  GLfloat *textureData = new GLfloat[4*nrOfVertices];

  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
  glReadPixels(0, 0, textureSize, textureSize,GL_RGBA, GL_FLOAT, &textureData[0]);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  ///---

  arcball_rotate();

  lightShader->use();
  volumeGenerator->render();
  lightShader->disable();

}

void WindowHandler::setupTextures()
{
    ///Example, if we need the position of the vertices
    //Get the position data for each vertex

	GLfloat *positionData = new GLfloat[volumeGenerator->getTetrahedMesh().GetVertexArraySize()];

	for(int i = 0; i < volumeGenerator->getTetrahedMesh().GetVertexArraySize(); i++)
	{
	    positionData[i] = volumeGenerator->getTetrahedMesh().GetVertexArray()[i];
	}


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

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void WindowHandler::init()
{

    lightShader = new Shader();
    lightShader->load("Shader/vertexPhongShader.glsl","Shader/fragmentPhongShader.glsl");

    deformationShader = new Shader();
    deformationShader->load("Shader/vertexDeformationShader.glsl","Shader/fragmentDeformationShader.glsl");

	volumeGenerator = new VolumeGenerator();
	volumeGenerator->generateVolume();

	//For the deformation
    textureSize = volumeGenerator->getTetrahedMesh().GetVertexArraySize(); //a texture is optimal if 2^n large
    nrOfVertices = textureSize*textureSize;



    //Setup textures used for the deformation shader
	setupTextures();




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
    gluPerspective( 45.0f, aspect_ratio, 0.1f, 1000.0f );
    gluLookAt(
        eye.x, eye.y, eye.z,   //eye
        center.x, center.y, center.z,   //lookat
        up.x, up.y, up.z );  //up vector
    arcball_setzoom( SPHERE_RADIUS, eye, up );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;

}
void WindowHandler::idle()
{

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
    if( buttonPressed == GLUT_LEFT_BUTTON || buttonPressed == GLUT_RIGHT_BUTTON || buttonPressed == GLUT_MIDDLE_BUTTON) {

		if( buttonPressed == GLUT_LEFT_BUTTON) {
			int invert_y = (windowHeight - y) - 1;
			arcball_move(x,invert_y);
		}

		glutPostRedisplay();

	}
}


void WindowHandler::processNormalKeys(unsigned char key, int x, int y) {

	if(key == 110) {

		volumeGenerator->changeNormalRenderMode();
	}

	if(key == 101) {

		volumeGenerator->changeEdgeRenderMode();
	}

	if(key == 116) {

		volumeGenerator->changeTriangleRenderMode();
	}

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



