#include <iostream>
#include <math.h>
#include "ScalarField.h"

template <class T>
ScalarField<T>::ScalarField(int xSize, 
                         int ySize, 
                         int xPadding, 
                         int yPadding, 
                         float xDisplacement, 
                         float yDisplacement,
                         float dx):xSize(xSize),ySize(ySize),xPadding(xPadding),yPadding(yPadding),xDisplacement(xDisplacement),yDisplacement(yDisplacement),dx(dx),count((xPadding*2+xSize)*(yPadding*2+ySize))
                        ,rowLength((xSize+xPadding*2))
                        ,invDx(1/dx){
    direction = NONE;
    //cout << "count: " << count;
    scalars = new T[count] ;
    new_scalars = new T[count];
                            
    //Fyll
    fill(scalars, scalars+(count),T());
    //Kopiera
    copy(scalars, scalars+count, new_scalars);
}

//Public
template<class T>
T ScalarField<T>::valueAtGridIndex(const int &x,const int &y){
    int i = indexAtGridIndex(x, y);
    //cout << "i :" << i << " x: " << x << " y: " << y << "count: " << count << endl;
    return scalars[i];
}


template<class T>
T ScalarField<T>::valueAtWorldCoord(float x,float y){
    
    int i, j;
    upperLeftCellIndexFromWorld(x, y, i, j);
    //static float xDis = xDisplacement*dx;
    //static float yDis = xDisplacement*dx;

    //i och j index uttryckta i världskoordinater...
    float iInWorld = (i+xDisplacement)*dx; 
    float jInWorld = (j+yDisplacement)*dx;
    //...och hitta q och k
    float q = (x-iInWorld)*invDx;
    float k = (y-jInWorld)*invDx;

    int index;
    
    
    index = indexAtGridIndex(i, j-1);
    T p1 = catmullRom(scalars[index-1], 
                          scalars[index],
                          scalars[index+1],
                          scalars[index+2], q);
    
    index = indexAtGridIndex(i, j);
    T p2 = catmullRom(scalars[index-1], 
                          scalars[index],
                          scalars[index+1],
                          scalars[index+2], q);
    
    index = indexAtGridIndex(i, j+1);
    T p3 = catmullRom(scalars[index-1], 
                          scalars[index],
                          scalars[index+1],
                          scalars[index+2], q);
    
    index = indexAtGridIndex(i, j+2);
    T p4 = catmullRom(scalars[index-1], 
                          scalars[index],
                          scalars[index+1],
                          scalars[index+2], q);
    
    T p = catmullRom(p1, p2, p3, p4, k);
    
    return p;
}

//Hitta array index för motsvarande gridindex
template<class T>
int ScalarField<T>::indexAtGridIndex(const int &x,const int &y){
    int i = (yPadding+y)*(rowLength)+(x+xPadding);
    return i;
}

//Hitta världskoordinaterna för cellindex
template<class T>
void ScalarField<T>::worldCoordinateAtCellIndex(const int &x, const int &y, float &world_x, float &world_y){
    world_x = (x+xDisplacement)*dx;
    world_y = (y+yDisplacement)*dx;
    //cout << "world_x: " << world_x << " world_y: " << world_y << endl;
}

//Hitta närmaste index från en världskoordinat
template<class T>
int ScalarField<T>::closestIndexFromWorld(const float &world_x, const float &world_y){
    float local_x = world_x-xDisplacement*dx;
    float local_y = world_y-yDisplacement*dx;

    int i = roundf(local_x*invDx);
    int j = roundf(local_y*invDx);
    
    
    return indexAtGridIndex(i, j); 
}

//Hitta närmaste index från en världskoordinat
template<class T>
int ScalarField<T>::closestIndexFromWorld(float world_x, float world_y, bool saveModeOn){

    if (saveModeOn) {
        world_x = world_x < 1 ? 1 : ( world_x > (xSize-1)*dx-1 ? xSize*dx-1 : world_x );
        world_y = world_y < 1 ? 1 : ( world_y > (ySize-1)*dx-1 ? ySize*dx-1 : world_y );
    }
    
    return closestIndexFromWorld(world_x, world_y); 
}

template<class T>
T ScalarField<T>::maxValue(){
    if(count < 0) return NULL;
    
    T max = scalars[0];
    for (int i = 1; i < count; i++) {
        if (max > scalars[i]) {
            max = scalars[i];
        }
    }
    
    return max;
}

//Substract
template<class T>
void ScalarField<T>::addValueAtIndex(T valueToAdd, const int &index){
    scalars[index] += valueToAdd;
    new_scalars[index] += valueToAdd;
}

//Sätt värde för ett index
template<class T>
void ScalarField<T>::setValueAtIndex(T newValue, const int &index){
    scalars[index] = T(newValue);
    new_scalars[index] = T(newValue);
}

//Sätt ett värde som kommer användas vid nästa bufferbyte
template<class T>
void ScalarField<T>::setFutureValueAtIndex(T newValue, const int &index){
    new_scalars[index] = newValue;
}

template<class T>
void ScalarField<T>::swapBuffer(){
    T *temp = scalars;
    scalars = new_scalars;
    new_scalars = temp;
}

//Private
template<class T>
T ScalarField<T>::catmullRom(const T &p1 ,const T &p2, const T &p3, const T &p4, const float &q){
    float q2 = q*q;
    float q3 = q2*q;
    return p1*(-0.5f*q+q2-0.5f*q3)+p2*(1.0f-2.5f*q2+1.5f*q3)+p3*(0.5f*q+2*q2-1.5f*q3)+p4*(-0.5f*q2+0.5f*q3);
}

template<class T>
int ScalarField<T>::upperLeftIndexFromWorld(float x, float y){
    int i = 0, j = 0;    
    upperLeftCellIndexFromWorld(x, y, i, j);
    return indexAtGridIndex(i, j);
}

template<class T>
void ScalarField<T>::upperLeftCellIndexFromWorld(float x,float y,int &i, int &j){
    //Gör om till lokala koordinater
    x -= xDisplacement*dx; 
    y -= yDisplacement*dx;
    
    //cout << "vx: " << x << " vy: " << y << endl;
    
    //Hitta närmaste index upp till vänster
    i = floorf(x*invDx);
    j = floorf(y*invDx);
}

template<class T>
ScalarField<T>::~ScalarField(){
    delete [] scalars;
    delete [] new_scalars;
}

template<class T>
ostream& operator << (ostream& ostr, const ScalarField<T>& d){
    for (int i = 0; i < d.count; i++) {
        ostr << " ";
        ostr.width(4);
        ostr << d.scalars[i] << endl;
    }
    return ostr;
}








