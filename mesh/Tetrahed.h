#include <vector>
using namespace std;

class Tetrahed {
    public:
        Tetrahed();
        ~Tetrahed();
        unsigned int getVertexInd();
        unsigned int getEdgeInd();
        vector<unsigned int> getFaceInd();
        

        
        void setFaceInd1(unsigned int f);
		void setFaceInd2(unsigned int f);
		void setFaceInd3(unsigned int f);
		void setFaceInd4(unsigned int f);
       

    private:
       
        unsigned int faceInd1; //face
		unsigned int faceInd2; //face
		unsigned int faceInd3; //face
		unsigned int faceInd4; //face
      
       
};
