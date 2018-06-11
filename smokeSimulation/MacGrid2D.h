
#ifndef smokeSimulation_MacGrid2D_h
#define smokeSimulation_MacGrid2D_h

#include "ScalarField.h"
#include "ScalarField.cpp"
#include "Vector2D.h"
#include "Vector2D.cpp"


enum Type {
    SOLID = 0,
    FLUID = 1,
    AIR = 2};


class MacGrid2D {
    //Funktioner
public:
    MacGrid2D(int xSize,int ySize, float dx, float dt);
    ~MacGrid2D();

    void RK2(const float &pos_x, const float &pos_y, float &new_x_pos, float &new_y_pos);
    void advect();
    void externalForces();
    void fastProject();
    //Lös tryckekvationen Ax = rhs...
    double* LGS(double *rhs, double *b, int rhsSize, int iterations, double limit, double scale); //...returnera x
    
    
    //Hantera skalärfält
    void initFields();
    void extrapolate();
    
    //Particles;
    void initParticles();
    void moveParticles();

    //Current steps
    int getSteps();
    
    //Variabler
public:

    const int xSize,ySize;
    const float dx;
    int stepCounter;
    float dt;
    float rho;
    ScalarField<Type> types;
    SootField<float> temprature;
    ScalarField<float> u, v;
    ScalarField<Type> cell;
    SootField<float> soot;
    //ScalarField<float> water;

    const int particleCount;
    Vector2D<float> *particles;
    double *pressure;
};

#endif
