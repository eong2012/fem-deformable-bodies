#include "Shader.h"
#include <iostream>
using namespace std;

Shader::Shader()
{
    programObject = 0;
    vertexShader = 0;
    fragmentShader = 0;

    uniformLocationTime = -1;
}

Shader::~Shader()
{
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    glDeleteProgram(programObject);
}

// Print a error
void Shader::printError( const char *errtype, const char *errmsg )
{
    fprintf( stderr, "%s: %s\n", errtype, errmsg );
}

// Get the filelength of the shader file
int Shader::filelength( const char *filename )
{
    FILE *ifp;
    int length = 0;

    ifp = fopen( filename, "r" );
    fseek( ifp, 0, SEEK_END );
    length = (int)ftell( ifp );
    fclose( ifp );
    return length;
}

// Read the shader file to a string
char* Shader::read(const char *filename)
{
    FILE *file = fopen( filename, "r" );
    if ( file == NULL )
    {
        printError( "I/O error", "Cannot open shader file!" );
        return 0;
    }
    int bytesinfile = filelength( filename );
    char *buffer = (char*)malloc( bytesinfile+1 );
    int bytesread = fread( buffer, 1, bytesinfile, file );
    buffer[bytesread] = 0; // Terminate the string with 0
    fclose( file );

    return buffer;
}

// Load the vertex shader and the fragment shader
void Shader::load(const char *vertexfile, const char *fragmentfile)
{
    GLint vertexCompiled;
    GLint fragmentCompiled;
    GLint shadersLinked;

    char *vertexSrcCode = NULL;
    char *fragmentSrcCode = NULL;

    char str[4096]; // For error messages from the GLSL compiler and linker

    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    vertexSrcCode = read( vertexfile );
    fragmentSrcCode = read( fragmentfile );

    const char *vCode = vertexSrcCode;
    const char *fCode = fragmentSrcCode;

    glShaderSource(vertexShader, 1, &vCode, NULL);
    glShaderSource(fragmentShader, 1, &fCode, NULL);

    free(vertexSrcCode);
    free(fragmentSrcCode);

    /* Compile the vertex shader, get the compile status, and print the log if something goes wrong */
    glCompileShader(vertexShader);

    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &vertexCompiled);

    if (vertexCompiled == GL_FALSE)
    {
        glGetShaderInfoLog( vertexShader, sizeof(str), NULL, str );
        printError( "Vertex shader compile error", str );
    }

    /* Compile the fragment shader, get the compile status and print the log if something went wrong */
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &fragmentCompiled);

    if (fragmentCompiled == GL_FALSE)
    {
        glGetShaderInfoLog( vertexShader, sizeof(str), NULL, str );
        printError( "Fragment shader compile error", str );
    }

    programObject = glCreateProgram();

    glAttachShader(programObject, vertexShader);
    glAttachShader(programObject, fragmentShader);

    /* Link the shaders, get link status and print the log if something went wrong */
    glLinkProgram(programObject);

    glGetProgramiv(programObject, GL_LINK_STATUS, &shadersLinked);

    if (shadersLinked == GL_FALSE )
    {
        glGetProgramInfoLog( programObject, sizeof(str), NULL, str );
        printError( "Program object linking error", str );
    }


}
void Shader::use()
{
    glUseProgram(programObject);
}
void Shader::disable()
{
    glUseProgram(0);
}

// Set texture to shader
void Shader::sendUniformTexture(const char *name,int unitIndex)
{

    GLint uniformLocation = glGetUniformLocation(programObject,(const char*) name);

    if (uniformLocation == -1)
        printError("Binding error","Failed to locate uniform texture");

    glUniform1i(uniformLocation,unitIndex);

}

void Shader::sendUniformFloat(const char *name,float f)
{
    GLint uniformLocation = glGetUniformLocation(programObject,(const char*) name);

    if (uniformLocation == -1)
        printError("Binding error","Failed to locate uniform variable");

    glUniform1f( uniformLocation, f );
}




