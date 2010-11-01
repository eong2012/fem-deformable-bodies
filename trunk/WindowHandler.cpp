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


void WindowHandler::display()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    Render();

    glutSwapBuffers();
}

void WindowHandler::Render()
{
	glEnable(GL_CULL_FACE);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 // glMatrixMode(GL_MODELVIEW);
  //glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  arcball_rotate();
  shader->use();
 
   
  volumeGenerator->render();
  shader->disable();
   //glClear(GL_COLOR_BUFFER_BIT);

//	glFlush();

}


void WindowHandler::init()
{
	
    shader = new Shader();
    shader->load("Shader/vertexPhongShader.glsl","Shader/fragmentPhongShader.glsl");
	volumeGenerator = new VolumeGenerator();
	volumeGenerator->generateVolume();


   
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



