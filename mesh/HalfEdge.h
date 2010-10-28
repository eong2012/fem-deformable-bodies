
class HalfEdge {
    public:
        HalfEdge();
        ~HalfEdge();
        unsigned int getVertexInd();
        unsigned int getFaceInd();
        unsigned int getPairInd();
        unsigned int getNextInd();
        unsigned int getPrevInd();

        void setVertexInd(unsigned int v);
        void setFaceInd(unsigned int f);
        void setPairInd(unsigned int p);
        void setNextInd(unsigned int n);
        void setPrevInd(unsigned int p);

    private:
        unsigned int vertexInd; //vertex
        unsigned int faceInd; //face
        unsigned int pairInd; //edge
        unsigned int nextInd; //edge
        unsigned int prevInd; //edge


};
