#include <iostream>
#include "Vector2D.h"
#include "math.h"
#include <fstream>

using namespace std;

template <class T>
Vector2D<T>::Vector2D():x(0), y(0){

}

template <class T>
Vector2D<T>::Vector2D(T newx, T newy):x(newx),y(newy){
    
}

template <class T>
Vector2D<T> Vector2D<T>::getNormalized(){
    float invSize = 1/length();
    return Vector2D(x*invSize,y*invSize);
}

template <class T>
T Vector2D<T>::length(){
    return sqrt(x*x+y*y);
}

template <class T>
void Vector2D<T>::operator+= (const Vector2D<T> &v){
    x += v.x;
    y += v.y;
}

template <class T>
Vector2D<T> Vector2D<T>::operator+ (const Vector2D<T> &v){
    return Vector2D(x+v.x,y+v.y);
}

template <class T>
Vector2D<T> Vector2D<T>::operator* (const T &scale){
    return Vector2D(x*scale,y*scale);
}

template <class T>
ostream& operator << (std::ostream& ostr, const Vector2D<T>& d){
    
    ostr << " (";
    ostr.fill(' ');
    ostr.width(4);
    ostr << d.x;
    ostr << " , ";
    ostr.fill(' ');
    ostr.width(4);
    ostr << d.y << ") "; 
    return ostr;
}