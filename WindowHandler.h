#include "Shader/Shader.h"
#include <GL/glew.h>
#include <GL/glut.h>

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>



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


private:


    int windowWidth, windowHeight;
    int textureSize;

    Shader *shader;

};
#endif
