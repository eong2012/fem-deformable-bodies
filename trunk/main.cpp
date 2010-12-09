#include "WindowHandler.h"
#include "AntTweakBar.h"
#include "../usr/include/armadillo.h"

using namespace std;

WindowHandler *mainWindow = 0;


void display(void)                              { mainWindow->display(); 	                          }
void idle(void)                                 { mainWindow->idle();                               }
void reshape(int width, int height)             { mainWindow->reshape(width, height);               }
void mouse(int button, int state, int x, int y) { mainWindow->mouseButtonEvent(button,state,x,y);   }
void move(int x, int y)                         { mainWindow->mouseMoveEvent(x,y);                  }
void processNormalKeys(unsigned char key, int x, int y)	{ mainWindow->processNormalKeys(key, x, y);		}

int main(int argc, char* argv[])
{




	glutInit(&argc, argv);
	TwInit(TW_OPENGL, NULL);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

	mainWindow = new WindowHandler();

	GLenum err = glewInit();



    mainWindow->init();
    glClearColor(1., 1., 1., 0.);

	glDisable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
    glutReshapeFunc(reshape);
	/*glutKeyboardFunc(processNormalKeys);
    glutMouseFunc(mouse);
    glutMotionFunc(move);*/

	 // Set GLUT event callbacks
    // - Directly redirect GLUT mouse button events to AntTweakBar
    glutMouseFunc((GLUTmousebuttonfun)TwEventMouseButtonGLUT);
    // - Directly redirect GLUT mouse motion events to AntTweakBar
    glutMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
    // - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
    glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
    // - Directly redirect GLUT key events to AntTweakBar
    glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
    // - Directly redirect GLUT special key events to AntTweakBar
    glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
    // - Send 'glutGetModifers' function pointer to AntTweakBar;
    //   required because the GLUT key event functions do not report key modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);





    //TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with GLUT and OpenGL.' "); // Message added to the help bar.
    //TwDefine(" TweakBar size='60 60' color='96 216 224' "); // change default tweak bar size and color



	glutIdleFunc(idle);
	glutMainLoop();
	TwTerminate();
	return 0;
}
