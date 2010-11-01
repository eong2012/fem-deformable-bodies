#include "Shader/Shader.h"
#include "mesh/HalfFaceMesh.h"

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
    void Render();

    //Functions for controlling the view
    void mouseButtonEvent(int button, int state, int x, int y);
    void mouseMoveEvent(int x, int y);


private:


    int windowWidth, windowHeight;
    int textureSize;

    HalfFaceMesh *halfFaceMesh;

    Shader *shader;

    //Arcball stuff
    vec eye;
    vec center;
    vec up;
    float SPHERE_RADIUS;
    float PI;
    int buttonPressed;


};
#endif
