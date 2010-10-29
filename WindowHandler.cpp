#include "WindowHandler.h"

WindowHandler::WindowHandler(void)
{

    windowWidth = 600;
    windowHeight = 600;
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("");

}
WindowHandler::~WindowHandler(void)
{

}



void WindowHandler::display()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);

    Render();
    glutSwapBuffers();
}

void WindowHandler::Render()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 // glMatrixMode(GL_MODELVIEW);
  //glLoadIdentity();

  halfFaceMesh->Render();
   //glClear(GL_COLOR_BUFFER_BIT);

//	glFlush();



}


void WindowHandler::init()
{

  //  shader = new Shader();
  //  shader->load("Shader/vertex.glsl","Shader/fragment.glsl");

    halfFaceMesh = new HalfFaceMesh();
    halfFaceMesh->AddTetrahedron();
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
        0.0f, 0.0f, 2.0f,   //eye
        0.0f, 0.0f, 0.0f,   //lookat
        0.0f, 1.0, 0.0f );  //up vector

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;

}
void WindowHandler::idle()
{

    glutPostRedisplay();
}



