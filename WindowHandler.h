#include "Shader/Shader.h"
#include <GL/glew.h>
#include <GL/glut.h>

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

#include "GLstuff/GLObject.h"
#include "Geometry/Mesh.h"
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
    bool Load(Mesh *mesh); // false return on error
    bool ReadHeader(std::istream &is);
    bool ReadData(std::istream &is);

    static Vector3<unsigned int> ReadTri(std::istream &is);

    template <class MeshType>
    MeshType * WindowHandler::AddMesh()
    {
    /*  String filename = path.AfterLast('/');
      if (filename == path) // If we're on Windows
        filename = path.AfterLast('\\');
      String suffix = path.AfterLast('.');

      if (suffix == _T("obj")) {*/
        // Create new mesh
        MeshType * mesh = new MeshType();
        mesh->SetName("sphereSmall");

        // Load mesh and add to geometry list
        std::ifstream infile;
        ObjIO objIO;
        infile.open("Mesh/sphereSmall.obj/");
        objIO.Load(mesh, infile);
        mesh->Initialize();

        // Add mesh to scene
        AddUniqueObject(mesh);

        return mesh;
      /*}
      else
        std::cerr << "Error: File type not supported" << std::endl;*/

      return NULL;
    }

private:


    int windowWidth, windowHeight;
    int textureSize;

    Shader *shader;

    struct LoadData{
    std::vector<Vector3<float> > verts;
    std::vector<Vector3<unsigned int> > tris;
  } loadData;

};
#endif
