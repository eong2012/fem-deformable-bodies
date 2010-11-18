
#define GLUT_DISABLE_ATEXIT_HACK
#define GLEW_STATIC
#ifndef SHADER_H
#define SHADER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>

class Shader
{
public:
	Shader();
	~Shader();

	void load(const char *vertexfile, const char *fragmentfile);
    void use();
    void disable();

	void sendUniformTexture(const char *name,int unitIndex);
	void sendUniformFloat(const char *name, float f);
	void initDataTexture(GLuint *texID, float* speeds, int maxParticles, int activeNr);

private:

	GLuint programObject;
	GLuint vertexShader;
	GLuint fragmentShader;

    GLint uniformLocationTime;



	void printError( const char *errtype, const char *errmsg );
	int filelength( const char *filename );
	char* read(const char *filename);


};
#endif
