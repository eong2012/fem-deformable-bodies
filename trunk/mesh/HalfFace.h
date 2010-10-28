class HalfFace {
    public:
        HalfFace();
        ~HalfFace();
        unsigned int getVertexInd();
        unsigned int getEdgeInd();
        unsigned int getFaceInd();
        unsigned int getMateInd();
        unsigned int getNextInd();
        unsigned int getPrevInd();
        unsigned int getRadialInd();

        void setVertexInd(unsigned int v);
        void setEdgeInd(unsigned int e);
        void setFaceInd(unsigned int f);
        void setMateInd(unsigned int m);
        void setNextInd(unsigned int n);
        void setPrevInd(unsigned int p);
        void setRadialInd(unsigned int r);

    private:
        unsigned int vertexInd; //vertex
        unsigned int edgeInd; //edge
        unsigned int faceInd; //face
        unsigned int mateInd; //edge
        unsigned int nextInd; //edge
        unsigned int prevInd; //edge
        unsigned int radialInd; //edge on the opposite face
};
