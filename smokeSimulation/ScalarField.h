#ifndef smokeSimulation_ScalarField_h
#define smokeSimulation_ScalarField_h


enum Direction {
    NONE = 0,
    VERTICAL = 1,
    HORIZONTAL = 2,};

enum Boundary {
    LEFT = 1,
    RIGHT = 2,
    TOP = 3,
    BOTTOM = 4,
    };

using namespace std;
template <class T>
class ScalarField {
    
    //Funktioner
public:
    ScalarField(int xSize, int ySize, int yPadding, int xPadding, float xDisplacement, float yDisplacement, float dx);
    void initScalars();
    ~ScalarField();
    
    //Index
    T valueAtGridIndex(const int &x,const int &y);
    T valueAtWorldCoord(float x, float y);
    T valueAtWorldCoord(float x, float y, bool saveModeOn);
    int indexAtGridIndex(const int &x,const int &y);
    void worldCoordinateAtCellIndex(const int &x, const int &y, float &world_x, float &world_y);
    int closestIndexFromWorld(const float &world_x, const float &world_y);
    int closestIndexFromWorld(float world_x, float world_y, bool saveModeOn);

    //Värden
    T maxValue();
    void addValueAtIndex(T valueToAdd, const int &index);
    void setValueAtIndex(T newValue, const int &index);
    void setFutureValueAtIndex(T newValue, const int &index); //Till buffer..
    void swapBuffer();
    
    
    //friend ostream& operator << (ostream& ostr, const ScalarField<T>& d);
protected:
    T catmullRom(const T &p1 ,const T &p2, const T &p3, const T &p4, const float &x);
    int upperLeftIndexFromWorld(float x, float y);
    void upperLeftCellIndexFromWorld(float x,float y,int &i, int &j);
    //Variabler
public:
    const int xSize, ySize, yPadding, xPadding, count;
    const float xDisplacement, yDisplacement;
    const float dx;
    const float invDx;
    Direction direction;
private:
    const int rowLength;
protected:
    T *scalars;
    T *new_scalars;
};

template <class T>
class SootField : public ScalarField<T> {
    public:
    SootField<T>(int xSize, int ySize, int yPadding, int xPadding, float xDisplacement, float yDisplacement, float dx):ScalarField<T>(xSize,ySize,yPadding,xPadding,xDisplacement,yDisplacement,dx){
        min = 0.01f;
        max = 1.0f;
    };
    
    T min, max;

    //Sätt värde för ett index
    void setValueAtIndex(float newValue, const int &index){
        newValue = newValue > max ? max : (newValue < min ? min : newValue);
        ScalarField<T>::setValueAtIndex(newValue, index);
    };
    
    void addValueAtIndex(T valueToAdd, const int &index){
        valueToAdd = valueToAdd > max ? max : (valueToAdd < min ? min : valueToAdd);
        ScalarField<T>::addValueAtIndex(valueToAdd, index);
    };
    
    float valueAtWorldCoord(float x, float y){ return ScalarField<T>::valueAtWorldCoord(x, y); }
    
    float valueAtWorldCoord(float x, float y, bool saveModeOn){
        float value = ScalarField<T>::valueAtWorldCoord(x, y, saveModeOn);
        return value > max ? max : (value < min ? min : value);
    }
 
    
    //Sätt framtida värde
    void setFutureValueAtIndex(float newValue, const int &index){
        newValue = newValue > max ? max : (newValue < min ? min : newValue);
        ScalarField<T>::setFutureValueAtIndex(newValue, index);
    };
};



#endif
