
#include "Shader/Shader.h"
#include "mesh/VolumeGenerator.h"
#include "Solver/Solver.h"
#include <GL/glew.h>
#include <GL/glut.h>

#include <iostream>

#include <stdlib.h>
#include <stdio.h>


#include "arcball.h"

#include "../AntTweakBar/include/AntTweakBar.h"

using namespace std;

#ifndef WindowHandler_H
#define WindowHandler_H

class WindowHandler
{
public:
    WindowHandler(void);
    ~WindowHandler(void);

    void display();
    void idle();
    void init();
    void reshape(int width, int height);

    void RenderFirstPass();
    void RenderSecondPass();
    void setupTextures();

    void drawQuad();
	void drawForceArrow();

	void processNormalKeys(unsigned char key, int x, int y);

    //Functions for controlling the view
    void mouseButtonEvent(int button, int state, int x, int y);
    void mouseMoveEvent(int x, int y);

	void MultiplyQuaternions(const float *q1, const float *q2, float *qout);
	void ConvertQuaternionToMatrix(const float *quat, float *mat);
	void SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat);
	void keyHandler();
	friend void TW_CALL normalCB(void *clientdata);
	friend void TW_CALL edgeCB(void *clientdata);
	friend void TW_CALL triangleCB(void *clientdata);


private:

    int windowWidth, windowHeight;

	Solver* solver;
    VolumeGenerator *volumeGenerator;
    //Stuff for the second pass
    Shader *lightShader;

    //Stuff for the first pass
    GLuint positionTexID; //Add necessary textures
    GLuint fbo;
    Shader *deformationShader;
    int textureSize;
    int nrOfVertices;
    //Arcball stuff
    vec eye;
    vec center;
    vec up;
    float SPHERE_RADIUS;
    float PI;
    int buttonPressed;
	arma::Mat<double> Fxt;

	TwBar *bar;
	float g_ForceDirection[3];

	float g_Rotation[4];



	float g_RotateStart[4];

	float g_Force;
	int g_NormalMode;
	int g_EdgeMode;
	int g_TriangleMode;


	float g_Mass;
	float g_Dampening;

	float g_fractureThresh;

};
#endif
