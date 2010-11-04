#include "WindowHandler.h"

using namespace std;



WindowHandler *mainWindow = 0;


void display(void)                              { mainWindow->display();                            }
void idle(void)                                 { mainWindow->idle();                               }
void reshape(int width, int height)             { mainWindow->reshape(width, height);               }
void mouse(int button, int state, int x, int y) { mainWindow->mouseButtonEvent(button,state,x,y);   }
void move(int x, int y)                         { mainWindow->mouseMoveEvent(x,y);                  }
void processNormalKeys(unsigned char key, int x, int y)	{ mainWindow->processNormalKeys(key, x, y);		}

int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

	mainWindow = new WindowHandler();

	GLenum err = glewInit();

    mainWindow->init();
    glClearColor(1., 1., 1., 0.);

	glDisable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
    glutReshapeFunc(reshape);
	glutKeyboardFunc(processNormalKeys);
    glutMouseFunc(mouse);
    glutMotionFunc(move);
	glutIdleFunc(idle);
	glutMainLoop();


  return 0;


	return 0;
}
