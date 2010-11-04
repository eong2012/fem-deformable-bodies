#include "Shader/Shader.h"
#include "mesh/VolumeGenerator.h"

#include <GL/glew.h>
#include <GL/glut.h>

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

#include "arcball.h"

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

	void processNormalKeys(unsigned char key, int x, int y);

    //Functions for controlling the view
    void mouseButtonEvent(int button, int state, int x, int y);
    void mouseMoveEvent(int x, int y);


private:

    int windowWidth, windowHeight;


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


};
#endif