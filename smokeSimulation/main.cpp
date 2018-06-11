#include <iostream>

//#include "Vector2D.h"
//#include "ScalarField.h"
#include "ScalarField.h"

#include "MacGrid2D.h"
//#include "MacGrid2D.cpp"

//#include "Vector2D.h"
//#include "Vector2D.cpp"

#include <math.h>
#if defined(__APPLE__) || defined(__MACOSX)
#include "glfw.h" // Include OpenGL Framework library
#else
#define GLFW_DLL
#define GLFW_BUILD_DLL
#include <GL\glfw.h>
#endif


// Frame counter and window settings variables
int frame      = 0, width     = 600, height      = 600;
int redBits    = 8, greenBits = 8,   blueBits    = 8;
int alphaBits  = 8, depthBits = 0,   stencilBits = 0;

//Mouse Positions
int xpos = 1, ypos = 1;
int newxpos = 1, newypos = 1;

//Simulationflags
bool simulationRunning = true;
bool showVectorField = true;
bool showParticles = false;

double t0 = 0.0;
int frames = 0;
char titlestring[200];

using namespace std;
void initGL(int width, int height);
void showFPS(MacGrid2D &grid);
template <class T>
void drawScalarField(ScalarField<T> &field);

template <class T>
void drawGrid(ScalarField<T> &field);
template <class T>
void drawBorder(ScalarField<T> &field);
void drawVectorField(MacGrid2D &grid);
void drawParticles(MacGrid2D &grid);

template <class T>
void drawSoot(ScalarField<T> &field);
void drawSolid(ScalarField<Type> &field);
void showFPS(MacGrid2D &grid) {
    
    double t, fps;
    
    // Get current time
    t = glfwGetTime();  // Gets number of seconds since glfwInit()
    // If one second has passed, or if this is the very first frame
    if( (t-t0) > 1.0 || frames == 0 )
    {
        fps = (double)frames / (t-t0);
        sprintf(titlestring, "2D Smoke Simulation (%.1f FPS), dt = %.2f", fps,grid.dt);
        glfwSetWindowTitle(titlestring);
        t0 = t;
        frames = 0;
    }
    frames ++;
}

void initGL(int width, int height){
	// ----- Window and Projection Settings -----
    
	// Set the window title
	//glfwSetWindowTitle("GLFW Basecode");
    
	// Setup our viewport to be the entire size of the window
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);
    
	// Change to the projection matrix, reset the matrix and set up orthagonal projection (i.e. 2D)
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, width, height, 0, 0, 1); // Paramters: left, right, bottom, top, near, far
    
	// ----- OpenGL settings -----
    
	glfwSwapInterval(1); 		// Lock to vertical sync of monitor (normally 60Hz, so 60fps)
    
	glEnable(GL_SMOOTH);		// Enable (gouraud) shading
    
	glDisable(GL_DEPTH_TEST); 	// Disable depth testing
    
	glEnable(GL_BLEND);		// Enable blending (used for alpha) and blending function to use
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
	glLineWidth(5.0f);		// Set a 'chunky' line width
    
	glEnable(GL_LINE_SMOOTH);	// Enable anti-aliasing on lines
    
	glPointSize(5.0f);		// Set a 'chunky' point size
    
	glEnable(GL_POINT_SMOOTH);	// Enable anti-aliasing on points
}

template <class T>
void drawScalarField(ScalarField<T> &field){
    
    if (field.direction == HORIZONTAL) {
        glColor4f(0.5f, 0.5f, 1.0f, 1.0f);
    }else if(field.direction == VERTICAL){
        glColor4f(1.0f, 0.5f, 0.5f, 1.0f);

    }
    for (int x = -field.xPadding; x < field.xSize+field.xPadding; x++) {
        for (int y = -field.yPadding; y < field.ySize+field.yPadding; y++) {
            
            glBegin(GL_LINE_STRIP);

            glVertex2d(field.dx*(x+field.xDisplacement), (y+field.yDisplacement)*field.dx);
            if (field.direction == VERTICAL) {
                glVertex2d(field.dx*(x+field.xDisplacement), (y+field.yDisplacement)*field.dx+field.valueAtGridIndex(x, y));
            }else if(field.direction == HORIZONTAL){
                glVertex2d(field.dx*(x+field.xDisplacement)+field.valueAtGridIndex(x, y), (y+field.yDisplacement)*field.dx);
            }
    
            
            
            glEnd();
        }
    }
}

void drawVectorField(MacGrid2D &grid){
    
    glPopMatrix();
    glLineWidth(2.0f);
    glColor4f(1.0f, 1.0f, 1.0f,0.8f);
    
    for (int x = 0; x < grid.xSize; x++) {
        for (int y = 0; y < grid.ySize; y++) {
            
            float world_x = x*grid.dx+grid.dx*0.5f;
            float world_y = y*grid.dx+grid.dx*0.5f;

            float u = grid.u.valueAtWorldCoord(world_x, world_y);
            float v = grid.v.valueAtWorldCoord(world_x, world_y);
            //float size = 1/sqrt(u*u+v*v);
            
            
            
            glBegin(GL_LINE_STRIP);
            glVertex2d(world_x, world_y);
            glVertex2d(world_x+u*3, world_y+v*3);
            glEnd();
            /*            
            glPushMatrix();
            glTranslatef(world_x-0.5,world_y-0.5,0);
            
            static float toDegrees = 360.0f/(2*3.14);
            
            if (u != 0) {
                glRotatef(atanf(v/u)*toDegrees, 0.0, 0.0, 1.0);
            }
            
            glScalef(0.3*size, size, 0);

            glBegin(GL_LINE_STRIP);
            glVertex2d(0, 0);
            glVertex2d(0, 1);
            glVertex2d(-0.3, 1);
            glVertex2d(0.5, 1.5);
            glVertex2d(1.3, 1);
            glVertex2d(1, 1);
            glVertex2d(1, 0);
            glVertex2d(0, 0);
            glEnd();
        
            glPopMatrix();
             */
        }
    } 
}

void drawParticles(MacGrid2D &grid){
    glPopMatrix();
    glLineWidth(3.0f);
    glColor3f(1.0f, 0.1f, 0.1f);
    for (int i = 0; i < grid.particleCount; i++) {
        glBegin(GL_LINE_STRIP);
        glVertex2d(grid.particles[i].x, grid.particles[i].y);
        glVertex2d(grid.particles[i].x+2,grid.particles[i].y);
        glEnd();
    }
}

template <class T>
void drawBorder(ScalarField<T> &field){
    
    glLineWidth(1.0f);
    glColor4f(1.0f, 1.0f, 1.0f,0.7f);
    
    //Vågräta streck
    for (int y = 0; y < field.ySize+1; y += field.ySize) {
        glBegin(GL_LINE_STRIP);
        glVertex2d(0, y*field.dx);
        glVertex2d((field.xSize-1)*field.dx, y*field.dx);
        glEnd();
    }
    
    //lodräta streck
    for (int x = 0; x < field.xSize; x += field.xSize-1) {
        glBegin(GL_LINE_STRIP);
        glVertex2d(x*field.dx, 0);
        glVertex2d(x*field.dx, field.dx*field.ySize);
        glEnd(); 
    }
}

template <class T>
void drawGrid(ScalarField<T> &field){
    
    glLineWidth(1.0f);
    glColor4f(1.0f, 1.0f, 1.0f,0.7f);
    
    //Vågräta streck
    for (int y = 0; y < field.ySize+1; y++) {
        glBegin(GL_LINE_STRIP);
        glVertex2d(0, y*field.dx);
        glVertex2d((field.xSize-1)*field.dx, y*field.dx);
        glEnd();
    }
    
    //lodräta streck
    for (int x = 0; x < field.xSize; x++) {
        glBegin(GL_LINE_STRIP);
        glVertex2d(x*field.dx, 0);
        glVertex2d(x*field.dx, field.dx*field.ySize);
        glEnd();    
    }
}

template <class T>
void drawSoot(ScalarField<T> &field){
    float dx2 = 0.5f*field.dx;

    for (int x = 0; x < field.xSize; x++) {
        for (int y = 0; y < field.ySize; y++) {
            float sootValue = field.valueAtGridIndex(x, y);
            glColor4f(sootValue,sootValue, sootValue, sootValue);
            
            if(sootValue > 0){
                
            
                //glColor4f(1.0, 1.0, 1.0, 1.0);
                float world_x,world_y;
                field.worldCoordinateAtCellIndex(x, y, world_x, world_y);
            
                glBegin(GL_QUADS);
                glVertex2f(world_x-dx2, world_y-dx2);
                glVertex2f(world_x+dx2, world_y-dx2);
                glVertex2f(world_x+dx2, world_y+dx2);
                glVertex2f(world_x-dx2, world_y+dx2);

                glEnd();
            }
        }
    } 
}
void drawSolid(ScalarField<Type> &field){
    float dx2 = 0.5f*field.dx;
    
    for (int x = 0; x < field.xSize; x++) {
        for (int y = 0; y < field.ySize; y++) {
            Type sootValue = field.valueAtGridIndex(x, y);
            glColor4f(0,1, 1, 1);
            
            if(sootValue == SOLID){
                
                
                //glColor4f(1.0, 1.0, 1.0, 1.0);
                float world_x,world_y;
                field.worldCoordinateAtCellIndex(x, y, world_x, world_y);
                
                glBegin(GL_QUADS);
                glVertex2f(world_x-dx2, world_y-dx2);
                glVertex2f(world_x+dx2, world_y-dx2);
                glVertex2f(world_x+dx2, world_y+dx2);
                glVertex2f(world_x-dx2, world_y+dx2);
                
                glEnd();
            }
        }
    }
    
}

template <class T>
void drawTemprature(ScalarField<T> &field){
    float dx2 = 0.5f*field.dx;
    
    for (int x = 0; x < field.xSize; x++) {
        for (int y = 0; y < field.ySize; y++) {
            float sootValue = field.valueAtGridIndex(x, y);
            glColor4f(sootValue,0, 0, 1);
            
            if(sootValue > 0){
                
                
                //glColor4f(1.0, 1.0, 1.0, 1.0);
                float world_x,world_y;
                field.worldCoordinateAtCellIndex(x, y, world_x, world_y);
                
                glBegin(GL_QUADS);
                glVertex2f(world_x-dx2, world_y-dx2);
                glVertex2f(world_x+dx2, world_y-dx2);
                glVertex2f(world_x+dx2, world_y+dx2);
                glVertex2f(world_x-dx2, world_y+dx2);
                
                glEnd();
            }
        }
    } 
}

int main (int argc, const char * argv[]){
    
    // Flags
	bool running = true;

	// Initialise glfw
	glfwInit();

    float translationx = 0, translationy = 0;
    int xSize = 50;
    int ySize = 50;
    int brushSize  = 20;

    
    cout << argc << endl;
    
    if (argc == 4) {
        
        xSize = atoi(argv[1]);
        ySize = atoi(argv[2]);
        brushSize  = atoi(argv[3]);
        
        cout << "xSize: " << xSize << " ySize: " << ySize << " brushSize: " << brushSize << endl;
    }
    
    //ScalarField<float> field = ScalarField<float>(1,1,1,1,1,1,1);

    float dx = (width-2*translationx)/xSize;
    MacGrid2D macGrid = MacGrid2D(xSize, ySize, dx, 0.05f);

    //Vector2D vector = Vector2D(0,0);
    //cout << vector << endl;
    
    /*
    ScalarField u = ScalarField(4, 3, 1, 2, 0.5, 0.0, 50);

    int x = 1;
    int y = 1;
    u.valueAtGridIndex(x, y);*/
  /*  ScalarField u = ScalarField(4, 3, 1, 2, 0.0, 0.5, 100);
    u.direction = HORIZONTAL;
    ScalarField v = ScalarField(3, 4, 2, 1, 0.5, 0.0, 100);
    v.direction = VERTICAL;
*/

    //u.valueAtWorldCoord(22, 7);
    //int index = u.indexAtGridIndex(0, 0);
    //ScalarField solids = ScalarField(3,3,0,0,0.0,0.0,100);

   	// Create a window
	if(!glfwOpenWindow(width, height, redBits, greenBits, blueBits, alphaBits, depthBits, stencilBits, GLFW_WINDOW)){
		cout << "Failed to open window!" << endl;
		glfwTerminate();
		return 0;
    }
	// Call our initGL function to set up our OpenGL options
	initGL(width, height);
    
	while (running == true){
        
        showFPS(macGrid);
        
        // Increase our frame counter
        frame++;
        
        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT);
        // Reset the matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(translationx, translationy, 0.0f);
        //Börja rita här
        drawBorder(macGrid.u);
        
        if (showVectorField) {
            drawVectorField(macGrid);
            drawGrid(macGrid.u);
        }else {
            //drawTemprature(macGrid.temprature);
            drawSoot(macGrid.soot);
            drawSolid(macGrid.types);

        }
        
        if (showParticles) {
            drawParticles(macGrid);
        }
        
        if (glfwGetMouseButton(0)) { 
            //double oldxpos = xpos;
            //double oldypos = ypos;
            glfwGetMousePos(&xpos, &ypos);

            xpos -= translationx;
            ypos -= translationy;
            //glfwDisable(GLFW_MOUSE_CURSOR);
            
            
            
            float dx = macGrid.soot.dx;
            for (float i = -brushSize; i < brushSize; i++) {
                for (float j = -brushSize; j < brushSize; j++) {

                    macGrid.soot.addValueAtIndex(brushSize/(i*i+j*j), macGrid.soot.closestIndexFromWorld(xpos+dx*i, ypos+dx*j, true));
                    //macGrid.temprature.addValueAtIndex(brushSize/(i*i+j*j)*1.f, macGrid.temprature.closestIndexFromWorld(xpos+dx*i, ypos+dx*j, true));
                }
            }
            /*
            float r = brushSize;
            for (float x = -r; x < r; x++) {
                float yMax = sqrtf(r*r-(x*x));
                for ( float y = -yMax; y < yMax; y++) {
                    //printf("%f\n",r/sqrtf(x*x+y*y));
                    macGrid.soot.addValueAtIndex(r/(x*x+y*y), macGrid.soot.closestIndexFromWorld(xpos+dx*x, ypos+dx*y, true));
                }
            }*/

           // macGrid.v.setValueAtIndex(-10, macGrid.v.closestIndexFromWorld(xpos, ypos, true));
            

        }else if(glfwGetMouseButton(1)){
            float oldxpos = xpos;
            float oldypos = ypos;
            glfwGetMousePos(&xpos, &ypos);
            xpos -= translationx;
            ypos -= translationy;
            //glfwDisable(GLFW_MOUSE_CURSOR);
            
            //float vel_u = macGrid.u.valueAtWorldCoord(xpos, ypos);
            //float vel_v = macGrid.v.valueAtWorldCoord(xpos, ypos);
            float val = 20.0f*((xpos-oldxpos)/dx);
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true)-1);
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true)-1);
           
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true)-1);
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true));
            
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true)-1);
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true)+1);
            
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true));
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true)-1);
            
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true));
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true));
            
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true));
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true)+1);
            
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true)+1);
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true)-1);
            
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true)+1);
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true));
            
            macGrid.u.setValueAtIndex(val, macGrid.u.closestIndexFromWorld(xpos, ypos,true)+1);
            macGrid.v.setValueAtIndex(val, macGrid.v.closestIndexFromWorld(xpos, ypos, true)+1);
            
            
            
        }else{
            glfwEnable(GLFW_MOUSE_CURSOR);
        }
        
        if (glfwGetKey(GLFW_KEY_ENTER)) {
            showVectorField = true;
            //simulationRunning = true;
        }else if(glfwGetKey(GLFW_KEY_SPACE)){
            showVectorField = false;
            //simulationRunning = false;
        }
        
        if (glfwGetKey(GLFW_KEY_UP)) {
            showParticles = true;
        }
        
        if (glfwGetKey(GLFW_KEY_DOWN)) {
            showParticles = false;
        }
        
        if (simulationRunning) {
            /*float dt = 0.5;
            grid.constructVelocityField();*/
            macGrid.moveParticles();
            macGrid.extrapolate();
            macGrid.advect();
            macGrid.externalForces();
            macGrid.fastProject();
             
            //macGrid.moveParticles();
             
            //simulationRunning = false;
            //macGrid.u.setValueAtIndex(10, macGrid.u.indexAtGridIndex(40, 40));
            //cout << macGrid.getSteps() << endl;
            
            //macGrid.project();
            //macGrid.u.scalars[macGrid.u.indexAtGridIndex(6, 6)] = 10;
            //macGrid.u.swapBuffer();
            //macGrid.v.swapBuffer();
        }
        
        
        // ----- Stop Drawing Stuff! ------
        glfwSwapBuffers(); // Swap the buffers to display the scene (so we don't have to watch it being drawn!)
        
		//exit if ESC was pressed or window was closed
        running = !glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam(GLFW_OPENED);  
	}
    
	glfwTerminate();
    
    return 0;
}

