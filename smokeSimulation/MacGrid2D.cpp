
#include <iostream>
#include "MacGrid2D.h"
#include <math.h>

//template<class T>
MacGrid2D::MacGrid2D(int xSize,int ySize, float dx,float dt):xSize(xSize),ySize(ySize),dx(dx),temprature(xSize,ySize, 2,2,0.5f,0.5f,dx),u(xSize+1,ySize, 1,2,0.0f,0.5f,dx),v(xSize,ySize+1, 2,1,0.5f,0.0f,dx),particleCount(xSize*ySize*(4)),cell(xSize,ySize, 2,2,0.5f,0.5f,dx), soot(xSize*5,ySize*5, 3,3,0.5f,0.5f,dx*0.2f),types(xSize,ySize, 2,2,0.5f,0.5f,dx),dt(dt){

    //float world_x, world_y;
    //u.worldCoordinateAtCellIndex(0, 0, world_x, world_y);
    
    pressure = NULL;
    
    rho = 1;
    
    stepCounter = 0;
    
    //Initiera fälten
    initFields();
    
    //Extrapolera vid kanterna på skalärfälten
    extrapolate();
    
    //Initiera Partiklar
    initParticles();

    
}

void MacGrid2D::initParticles(){
    
    particles = new Vector2D<float>[particleCount];
    int counter = 0;
    //Skapa och placera partiklar
    for (int i = 0; i < xSize; i++) {
        for (int j = ySize/2; j < ySize; j++) {
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < 2; l++) {
                    particles[counter] = Vector2D<float>( (int)((i+0.5-0.25+k*0.5)*dx),(int)((j+0.5-0.25+l*0.5)*dx));
                    counter++;
                }
            }
        }
    }
}

void MacGrid2D::moveParticles(){
    for (int i = 0; i < particleCount; i++) {
        
        Vector2D<float> &particlePosition = particles[i];
        
        float u_vel, v_vel;
        u_vel = u.valueAtWorldCoord(particlePosition.x, particlePosition.y);
        v_vel = v.valueAtWorldCoord(particlePosition.x, particlePosition.y);

        particlePosition += Vector2D<float>(u_vel, v_vel)*(dt*0.5);
        u_vel = u.valueAtWorldCoord(particlePosition.x, particlePosition.y);
        v_vel = v.valueAtWorldCoord(particlePosition.x, particlePosition.y);
        
        particlePosition += Vector2D<float>(u_vel, v_vel)*dt;

    }
}

void MacGrid2D::initFields(){
    
    
    
    
    
    //Sätt celltyper
    for (int i = -cell.xPadding; i < cell.xSize+cell.xPadding; i++) {
        for (int j = -cell.yPadding; j < cell.ySize+cell.yPadding; j++) {
            
            if(i >= 0 && i < xSize && j >= 0 && j < ySize){
                cell.setValueAtIndex(FLUID, cell.indexAtGridIndex(i, j));
            }else{
                cell.setValueAtIndex(SOLID, cell.indexAtGridIndex(i, j));
            }
            
        }
    }
    
    u.direction = HORIZONTAL;
    v.direction = VERTICAL;
    
    //Sätt u värden
    for (int i = -u.xPadding; i < u.xSize+u.xPadding; i++) {
        for (int j = -u.yPadding; j < u.ySize+u.yPadding; j++) {
            
            if (cell.valueAtGridIndex(i, j) == FLUID && cell.valueAtGridIndex(i-1, j) == FLUID){
                u.setValueAtIndex(0,u.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
            }else{
                u.setValueAtIndex(0,u.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
            }
        }
    }
    
    //Sätt v värden
    for (int i = -v.xPadding; i < v.xSize+v.xPadding; i++) {
        for (int j = -v.yPadding; j < v.ySize+v.yPadding; j++) {
            //v.setValueAtIndex(0, v.indexAtGridIndex(i, j));
            //v.setValueAtIndex(((rand() % 1000000)*0.000001-0.5)*100, v.indexAtGridIndex(i, j));
            
            if (cell.valueAtGridIndex(i, j) == FLUID && cell.valueAtGridIndex(i, j-1) == FLUID){
                v.setValueAtIndex(0,v.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
            }else{
                v.setValueAtIndex(0,v.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
            }
        }
    }
    
    
    //Sätt soot värden
    for (int i = -soot.xPadding; i < soot.xSize+soot.xPadding; i++) {
        for (int j = -soot.yPadding+1; j < soot.ySize+soot.yPadding; j++) {

            if (cell.valueAtGridIndex(i, j) == FLUID && cell.valueAtGridIndex(i, j-1) == FLUID){
                soot.setValueAtIndex(0,soot.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
            }else{
                soot.setValueAtIndex(0,soot.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
            }
        }
    }
    
    //Sätt tempratur-värden
    for (int i = -temprature.xPadding; i < temprature.xSize+soot.xPadding; i++) {
        for (int j = -temprature.yPadding+1; j < temprature.ySize+soot.yPadding; j++) {
            temprature.setValueAtIndex(0,temprature.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
            temprature.setFutureValueAtIndex(0,temprature.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
        }
    }
    
    //Sätt typer
    //Sätt tempratur-värden
    for (int i = -types.xPadding; i < types.xSize+soot.xPadding; i++) {
        for (int j = -types.yPadding+1; j < types.ySize+soot.yPadding; j++) {
            
            
            if ((i > 5 && i < 20 && j > 5 && j < 20)) {
                types.setValueAtIndex(SOLID,temprature.indexAtGridIndex(i, j));
                types.setFutureValueAtIndex(SOLID,types.indexAtGridIndex(i, j));//((rand() % 1000000)*0.000001-0.5)*100 ;
            }else{
                types.setValueAtIndex(AIR,temprature.indexAtGridIndex(i, j));
                types.setFutureValueAtIndex(AIR,temprature.indexAtGridIndex(i, j));
            }
            
            

        }
    }


}


void MacGrid2D::extrapolate(){
    
    //Extrapolera de övre u värden
    for (int i = 1; i < xSize; i++) {
        for (int j = 0; j < 1; j++) {
            for (int k = 0; k < u.yPadding; k++) {
                float uVal = u.valueAtGridIndex(i, (u.ySize-1));
                u.setValueAtIndex(uVal,u.indexAtGridIndex(i, (u.yPadding+u.ySize-1-k)));// ((rand() % 1000000)*0.000001-0.5)*100;
            }
        }
    }
    
    //Extrapolera de undre u värden
    for (int i = 1; i < xSize; i++) {
        for (int j = 0; j < 1; j++) {
            for (int k = 0; k < u.yPadding; k++) {
                float uVal = u.valueAtGridIndex(i, 0);
                u.setValueAtIndex(uVal,u.indexAtGridIndex(i, (-1-k)));// ((rand() % 1000000)*0.000001-0.5)*100;
            }
        }
    }
    
    //Sätt v värden
    for (int i = 0; i < 1; i++) {
        for (int j = 1; j < ySize; j++) {
            for (int k = 0; k < v.xPadding; k++) {
                float vVal = v.valueAtGridIndex(v.xSize-1, j);
                v.setValueAtIndex(vVal,v.indexAtGridIndex((v.xPadding+v.xSize-1-k),j));// ((rand() % 1000000)*0.000001-0.5)*100;
            }
        }
    }
    
    //Sätt v värden
    for (int i = 0; i < 1; i++) {
        for (int j = 1; j < ySize; j++) {
            for (int k = 0; k < v.xPadding; k++) {
                float vVal = v.valueAtGridIndex(0, j);
                v.setValueAtIndex(vVal,v.indexAtGridIndex(-1-k,j));// ((rand() % 1000000)*0.000001-0.5)*100;
            }
        }
    }
    
    

}

void MacGrid2D::advect(){
    int x, y;
    
    
#if defined(__APPLE__) || defined(__MACOSX)
#define CHUNKSIZE 1
#pragma omp parallel private(x,y)
#pragma omp for schedule(static, CHUNKSIZE)
#endif

    
    //Advect soot
    for(x = 0; x < soot.xSize; x++){

        for (y = 0; y < soot.ySize; y++) {
            float world_x , world_y;

            //Hitta världskoordinaterna för cellindex
            soot.worldCoordinateAtCellIndex(x, y, world_x, world_y);
            
            //RK2
            float new_u, new_v;
            RK2(world_x, world_y, new_u, new_v);
            
            float val = soot.valueAtWorldCoord(world_x-new_u, world_y-new_v);
            soot.setFutureValueAtIndex(val*1.0f, soot.indexAtGridIndex(x, y));
        }
    }
    
    //Advect temprature
    for(x = 0; x < temprature.xSize; x++){
        
        for (y = 0; y < temprature.ySize; y++) {
            float world_x , world_y;
            
            //Hitta världskoordinaterna för cellindex
            temprature.worldCoordinateAtCellIndex(x, y, world_x, world_y);
            
            //RK2
            float new_u, new_v;
            RK2(world_x, world_y, new_u, new_v);
            
            float val = soot.valueAtWorldCoord(world_x-new_u, world_y-new_v);
            temprature.setFutureValueAtIndex(val*1.0f, temprature.indexAtGridIndex(x, y));
        }
    }
    


    //Advect U
    for (x = 1; x <  xSize; x++) {
        for (y = 0; y < ySize; y++) {
            //.. do stuff
            //Advect U

            float world_x , world_y;
                
                //Hitta världskoordinaterna för cellindex
                u.worldCoordinateAtCellIndex(x, y, world_x, world_y);
                
                //RK2
                float new_u, new_v;
                RK2(world_x, world_y, new_u, new_v);
                
                int index = u.indexAtGridIndex(x, y);
                //cout << "new_u " << new_u << endl;
                u.setFutureValueAtIndex(new_u, index); 
            
        }
    }
    
    //Advect V
    for (x = 0; x <  xSize; x++) {
        for (y = 1; y < ySize; y++) {
            //.. do stuff
            
            float world_x , world_y;
            
            //Hitta världskoordinaterna för cellindex
            v.worldCoordinateAtCellIndex(x, y, world_x, world_y);
            
            //RK2
            float new_u, new_v;
            RK2(world_x, world_y, new_u, new_v);
            
            int index = v.indexAtGridIndex(x, y);
            v.setFutureValueAtIndex(new_v, index);
        }
    }
    
    
    //Byt buffer - "aktivera" värden som är satta med setFutureValueAtIndex(index) 
    u.swapBuffer();
    v.swapBuffer();
    soot.swapBuffer();
    temprature.swapBuffer();
}

void MacGrid2D::externalForces(){
    
    //Advect U
    /*
     for (int x = 0; x <  xSize; x++) {
        for (int y = 1; y < ySize-1; y++) {
            //.. do stuff
            int i = v.indexAtGridIndex(x, y);
            v.addValueAtIndex(-9.82*dt, i);
        }
    }    
    */
    int x, y;
    //Advect V
    for (x = 0; x <  xSize; x++) {
        for (y = 1; y < ySize; y++) {
            //.. do stuff
            

            if (types.valueAtGridIndex(x, y) == AIR && types.valueAtGridIndex(x, y+1) == AIR) {
                float world_x , world_y;
                
                //Hitta världskoordinaterna för cellindex
                v.worldCoordinateAtCellIndex(x, y, world_x, world_y);
                float val = temprature.valueAtWorldCoord(world_x, world_y);
                
                int index = v.indexAtGridIndex(x, y);
                v.addValueAtIndex(-val*0.1f, index);
                
                
            }else{
                int index = v.indexAtGridIndex(x, y);
                v.setValueAtIndex(0, index);
            }
        }
    }
    
    //Advect U
    for (x = 0; x <  xSize; x++) {
        for (y = 1; y < ySize; y++) {
            //.. do stuff
            
            
            if (types.valueAtGridIndex(x-1, y) == AIR && types.valueAtGridIndex(x, y) == AIR) {
                
            }else{
                int index = u.indexAtGridIndex(x, y);
                u.setValueAtIndex(0, index);
            }
        }
    }
}

void MacGrid2D::fastProject(){
    
    int bCount  = xSize*ySize;
    double *b = new double[bCount];

    double scale = 1.0/dx;
    int index = 0;
    
    double mean = 0;
    for (int x = 0; x < xSize; x++) {
        for (int y = 0; y < ySize; y++) {
            b[index] = -scale*(u.valueAtGridIndex(x+1, y)-u.valueAtGridIndex(x, y)+v.valueAtGridIndex(x, y+1)-v.valueAtGridIndex(x, y));
            //b(index) = 0;
            mean += b[index];
            index++;
        }
    }

    mean *= 1.0/(bCount);
    //Ta bort medelvärdet
    for (int i = 0; i < bCount; i++) {
        b[i] -= mean;
    }
    
    
    
    //LaplacianGaussSeidel
    double limit = 1e-6;
    int iterations = 100;
    scale = dt/(rho * dx*dx);
    
    pressure = LGS(b, pressure, bCount, iterations, limit, scale);
    
    scale = dt/(rho*dx);
    
    //Applicera trycket på u rep. v
    index = 0;
    int i,j;
    for (i = 0; i < xSize; i++) {
        for (j = 0; j < ySize; j++) {
            
            //cout << scale*pressure[index] << endl;
            //Lägg på trycket
            //if (FLUID == types.valueAtGridIndex(i, j)) {
                float p = (float)(scale*pressure[index]);
                u.addValueAtIndex(-p, u.indexAtGridIndex(i, j));
                u.addValueAtIndex(p, u.indexAtGridIndex(i+1, j));
                v.addValueAtIndex(-p, v.indexAtGridIndex(i, j));
                v.addValueAtIndex(p, v.indexAtGridIndex(i, j+1));
            //}

            index++;
        }
    }
    
    
    //TODO!
    //Sätt nollor på kanterna!!!!! (soliderna)
    for (int i = -cell.xPadding; i < cell.xSize+cell.xPadding; i++) {
        for (int j = -cell.yPadding; j < cell.ySize+cell.yPadding; j++) {
            
            
            Type thisType = (Type)cell.valueAtGridIndex(i, j);
            Type leftType = (Type)cell.valueAtGridIndex(i-1, j);
            Type downType = (Type)cell.valueAtGridIndex(i, j-1);
            
            //Kommentar
            if ((thisType == SOLID && leftType == FLUID)
                ||(cell.valueAtGridIndex(i, j) == FLUID && leftType == SOLID)) {
                u.setValueAtIndex(0, u.indexAtGridIndex(i, j));
            }
            
            if ((thisType == SOLID && downType == FLUID)
                ||(thisType == FLUID && downType == SOLID)) {
                v.setValueAtIndex(0, v.indexAtGridIndex(i, j));
            }
        }
    }
    
    
    
    delete [] b;
}


double* MacGrid2D::LGS(double *b, double *guess, int rhsSize, int maxIterations, double limit, double scale){
    

    //Längden på A-matrisens diagonal
    int size = xSize*ySize;
    
    
    double *p;
    //Allokera p
    if (!guess) {
        p = new double[size];
        //Fyll p med 0-or
        fill(p, p+size, 0);
    }else{
        p = guess;
    }

    //double error = 0;
    //Spara felet efter iterationerna
    int iterations;
    //int nextCheck = 0; //Nästa gång som divergensen ska kollas
    scale = 1.0/scale;
    
    for (iterations = 0; iterations < maxIterations; iterations++) {
        int counter = 0;
        for (int i = 0; i < xSize; i++) {
            for (int j = 0; j < ySize; j++) {
                
                double diagElement = 0;
                double b_val = b[counter]*scale;
                
                double p1 =0 ,p2 = 0, p3 = 0,p4 = 0;
                
                Type bottom, left, right,up;
                bottom = types.valueAtGridIndex(i, j+1);
                left = types.valueAtGridIndex(i-1, j);
                right = types.valueAtGridIndex(i+1, j);
                up = types.valueAtGridIndex(i, j-1);
                
                //Bottom
                if (xSize+counter < size && bottom == AIR) {
                    p1 = p[counter+xSize];
                    diagElement += 1;
                }
                
                //Top
                if (counter-xSize >= 0 && up == AIR) {
                    p2 = p[counter-xSize];
                    diagElement += 1;
                }

                //Left
                if (counter-1 >= 0 && left == AIR) {
                    p3 = p[counter-1];
                    diagElement +=1;
                }
                
                //Right
                if (counter+1 < size && right == AIR) {
                    p4 = p[counter+1];
                    diagElement += 1;
                }
                
                if (diagElement == 0) {
                    p[counter] = 0;
                }else{
                    p[counter] = (p1+p2+p3+p4+b_val)/diagElement;

                }
                
                
                counter++;
            }
        }
        
        /*if (nextCheck == iterations) {
            //Räkna ut felet
            error = 0;
            counter = 0;
            for (int i = 0; i < xSize; i++) {
                for (int j = 0; j < ySize; j++) {
                    
                    double diagElement = 0;                
                    double p1 =0 ,p2 = 0, p3 = 0,p4 = 0;
                    
                    //Bottom
                    if (xSize+counter < size) {
                        p1 = p[counter+xSize];
                        diagElement += 1;
                    }
                    
                    //Top
                    if (counter-xSize >= 0 ) {
                        p2 = p[counter-xSize];
                        diagElement += 1;
                    }
                    
                    //Left
                    if (counter-1 >= 0 ) {
                        p3 = p[counter-1];
                        diagElement +=1;
                    }
                    
                    //Right
                    if (counter+1 < size) {
                        p4 = p[counter+1];
                        diagElement += 1;
                    }
                    
                    error += pow(((diagElement*p[counter]-p1-p2-p3-p4)/scale-b[counter]),2.0f);
                    counter++;
                    
                }
            } 
            
            //Lösningen är tillräckligt nära : avbryt
            if (sqrt(error) < limit) {
                cout << "error: "  << error << " at iterations: " << iterations << endl;  
                return p;
            }
            
            
            nextCheck += 10;
        }*/
        
    }
    
    //cout << "error: "  << error << " at iterations: " << maxIterations << endl; 
    
    
    return p;
}



void MacGrid2D::RK2(const float &x_pos, const float &y_pos , float &new_u_vel, float &new_v_vel){
    float u_comp = u.valueAtWorldCoord(x_pos, y_pos);
    float v_comp = v.valueAtWorldCoord(x_pos, y_pos);
    
    float temp_x = x_pos, temp_y = y_pos;
    
    //Ett halvt tiddsteg bakåt
    temp_x -= u_comp*0.5f*dt;
    temp_y -= v_comp*0.5f*dt;

    //Hitta nya vektorn i den nya punkten
    u_comp = u.valueAtWorldCoord(temp_x, temp_y);
    v_comp = v.valueAtWorldCoord(temp_x, temp_y);
    
    //Gå ett helt tidsteg bakåt
    temp_x -= u_comp*dt;
    temp_y -= v_comp*dt;
    
    //Hitta vektorn i den punkten
    new_u_vel = u.valueAtWorldCoord(temp_x, temp_y);
    new_v_vel = v.valueAtWorldCoord(temp_x, temp_y);
}

int MacGrid2D::getSteps(){
    return stepCounter++;
}

MacGrid2D::~MacGrid2D(){
    delete [] particles;   
    delete [] pressure;

}

