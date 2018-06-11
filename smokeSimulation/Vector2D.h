#ifndef fluidSimulation_Vector2D_h
#define fluidSimulation_Vector2D_h

template<class T>
class Vector2D {
public:
    Vector2D();
    Vector2D(T x,T y);
    Vector2D getNormalized();
    T length();

    void operator+= (const Vector2D &v);
    Vector2D operator+ (const Vector2D &v);
    Vector2D operator* (const T &scale);

    //friend std::ostream& operator << (std::ostream&, const Vector2D<T>& d);
    void swapVector(T x,T y);
    T x,y;
    Vector2D *vector;
};



#endif
