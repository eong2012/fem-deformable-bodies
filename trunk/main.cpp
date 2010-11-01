#include "WindowHandler.h"
#include <unistd.h>
using namespace std;

WindowHandler *mainWindow = 0;


void display(void)                              { mainWindow->display();                            }
void idle(void)                                 { mainWindow->idle();                               }
void reshape(int width, int height)             { mainWindow->reshape(width, height);               }
void mouse(int button, int state, int x, int y) { mainWindow->mouseButtonEvent(button,state,x,y);   }
void move(int x, int y)                         { mainWindow->mouseMoveEvent(x,y);                  }

int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

	mainWindow = new WindowHandler();

	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
	/* Problem: glewInit failed, something is seriously wrong. */
	fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}
	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
	if (!GLEW_ARB_vertex_buffer_object)
	{
	//cerr << "VBO not supported\n";
	exit(1);
	}

    mainWindow->init();
    glClearColor(1., 1., 1., 0.);

	glDisable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(move);
	glutIdleFunc(idle);
	glutMainLoop();

	return 0;
}
