#include<bits/stdc++.h>
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
using namespace std;

struct point{
    double x, y, z;
    point(double x = 0.0, double y = 0.0, double z = 0.0){
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

// Global variables
struct point pos;   // position of the eye
struct point l;     // look/forward direction
struct point r;     // right direction
struct point u;     // up direction

struct point sphereVertices[100][100];

const double RADIUS = 0.57735 ; // 1/sqrt(3)
const double HEIGHT = 1.41421 ; // sqrt(2)
const double STEP_COUNT = 16.0; // sqrt(2)

double octScale = 1.0;
double octTranslate = 0.0;
double sphereScale = 0.0;
double sphereTranslate = 1.0;
double cylinderTranslate = 0.5;

double rotateAngle = 0;
bool isAxes = false;

/* Initialize OpenGL Graphics */
void initGL() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   // Black and opaque
    glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
}

/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes() {
    glLineWidth(1);
    glBegin(GL_LINES);
        glColor3f(1,0,0);   // Red
        // X axis
        glVertex3f(0,0,0);
        glVertex3f(100,0,0);

        glVertex3f(0,0,0);
        glVertex3f(-100,0,0);

        glColor3f(0,1,0);   // Green
        // Y axis
        glVertex3f(0,0,0);
        glVertex3f(0,100,0);

        glVertex3f(0,0,0);
        glVertex3f(0,-100,0);

        glColor3f(0,0,1);   // Blue
        // Z axis
        glVertex3f(0,0,0);
        glVertex3f(0,0,100);

        glVertex3f(0,0,0);
        glVertex3f(0,0,-100);
    glEnd();
}

// generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
void drawSphereSegment(int subdivision)
{
    const float DEG2RAD = acos(-1) / 180.0f;

    float n1[3];        // normal of longitudinal plane rotating along Y-axis
    float n2[3];        // normal of latitudinal plane rotating along Z-axis
    float v[3];         // direction vector intersecting 2 planes, n1 x n2
    float a1;           // longitudinal angle along Y-axis
    float a2;           // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for(unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(unsigned int j = 0; j < pointsPerRow; ++j)
        {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float scale = (1 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]))*RADIUS;
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;

            point p(v[0],v[1],v[2]);
            // add a vertex into array
            sphereVertices[i][j]= p ;
        }
    }

    for(int i=0;i<pointsPerRow-1;i++){
        for(int j=0;j<pointsPerRow-1;j++){
            glBegin(GL_QUADS);
                glVertex3f(sphereVertices[i][j].x, sphereVertices[i][j].y, sphereVertices[i][j].z);
                glVertex3f(sphereVertices[i+1][j].x, sphereVertices[i+1][j].y, sphereVertices[i+1][j].z);
                glVertex3f(sphereVertices[i+1][j+1].x, sphereVertices[i+1][j+1].y, sphereVertices[i+1][j+1].z);
                glVertex3f(sphereVertices[i][j+1].x, sphereVertices[i][j+1].y, sphereVertices[i][j+1].z);
            glEnd();
        }
    }

    
}

void drawSphere(){
    glPushMatrix();
    
    //+X side
    glPushMatrix();
    glColor3f(1,0,0);
    glTranslated(sphereTranslate, 0, 0);
    glScaled(sphereScale, sphereScale, sphereScale);
    drawSphereSegment(3);
    glPopMatrix();

    //-Z side
    glPushMatrix();
    glColor3f(0,0,1);
    glRotated(90,0,1,0);
    glTranslated(sphereTranslate, 0, 0);
    glScaled(sphereScale, sphereScale, sphereScale);
    drawSphereSegment(3);
    glPopMatrix();

    //-X side
    glPushMatrix();
    glColor3f(1,0,0);
    glRotated(180,0,1,0);
    glTranslated(sphereTranslate, 0, 0);
    glScaled(sphereScale, sphereScale, sphereScale);
    drawSphereSegment(3);
    glPopMatrix();

    //+Z side
    glPushMatrix();
    glColor3f(0,0,1);
    glRotated(270,0,1,0);
    glTranslated(sphereTranslate, 0, 0);
    glScaled(sphereScale, sphereScale, sphereScale);
    drawSphereSegment(3);
    glPopMatrix();

    //+Y side
    glPushMatrix();
    glColor3f(0,1,0);
    glRotated(90,0,0,1);
    glTranslated(sphereTranslate, 0, 0);
    glScaled(sphereScale, sphereScale, sphereScale);
    drawSphereSegment(3);
    glPopMatrix();

    //-Y side
    glPushMatrix();
    glColor3f(0,1,0);
    glRotated(-90,0,0,1);
    glTranslated(sphereTranslate, 0, 0);
    glScaled(sphereScale, sphereScale, sphereScale);
    drawSphereSegment(3);
    glPopMatrix();
    
    glPopMatrix();
}

void drawTriangle(){
    glBegin(GL_TRIANGLES);           // Begin drawing the pyramid with 4 triangles
        glVertex3f( 1.0, 0.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 1.0);
    glEnd();
}

void drawOctahedron(){
    glPushMatrix();
    
    glPushMatrix();
    glColor3f(1, 0, 1);
    glTranslated(octTranslate, octTranslate,octTranslate);
    glScaled(octScale, octScale, octScale);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor3f(0, 1, 1);
    glRotated(90,0,1,0);
    glTranslated(octTranslate, octTranslate,octTranslate);
    glScaled(octScale, octScale, octScale);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor3f(1, 0, 1);
    glRotated(180,0,1,0);
    glTranslated(octTranslate, octTranslate,octTranslate);
    glScaled(octScale, octScale, octScale);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor3f(0, 1, 1);
    glRotated(270,0,1,0);
    glTranslated(octTranslate, octTranslate,octTranslate);
    glScaled(octScale, octScale, octScale);
    drawTriangle();
    glPopMatrix();


    glRotated(180,1,0,0);


    glPushMatrix();
    glColor3f(1, 0, 1);
    glTranslated(octTranslate, octTranslate,octTranslate);
    glScaled(octScale, octScale, octScale);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor3f(0, 1, 1);
    glRotated(90,0,1,0);
    glTranslated(octTranslate, octTranslate,octTranslate);
    glScaled(octScale, octScale, octScale);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor3f(1, 0, 1);
    glRotated(180,0,1,0);
    glTranslated(octTranslate, octTranslate,octTranslate);
    glScaled(octScale, octScale, octScale);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor3f(0, 1, 1);
    glRotated(270,0,1,0);
    glTranslated(octTranslate, octTranslate,octTranslate);
    glScaled(octScale, octScale, octScale);
    drawTriangle();
    glPopMatrix();

    glPopMatrix();
}

void drawCylinderSegment(double height, double radius, int slices) {
    double tempx = radius, tempy = 0;
    double currx, curry;
    glBegin(GL_QUADS);
        for (int i = -slices/2; i < slices/2; i++) {
            double theta = i * 0.4 * M_PI / slices;
            currx = radius * cos(theta);
            curry = radius * sin(theta);

            glVertex3f(currx, curry, height/2);
            glVertex3f(currx, curry, -height/2);

            glVertex3f(tempx, tempy, -height/2);
            glVertex3f(tempx, tempy, height/2);

            tempx = currx;
            tempy = curry;
        }
    glEnd();
}

void drawCylinder(){
    glPushMatrix();
    glColor3f(1,1,0);

    glPushMatrix();
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(90,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(180,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(270,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();
    
    glRotated(90, 0, 0, 1);

    glPushMatrix();
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(90,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(180,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(270,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glRotated(90, 0, 0, 1);

    glPushMatrix();
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(90,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(180,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(270,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();


    glRotated(90, 0, 0, 1);

    glPushMatrix();
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(90,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(180,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPushMatrix();
    glRotated(270,1,0,0);
    glTranslated(cylinderTranslate, 0, cylinderTranslate);
    glRotated(-45,0,1,0);
    glScaled(sphereScale, sphereScale, octScale);
    drawCylinderSegment(HEIGHT, RADIUS, 100);
    glPopMatrix();

    glPopMatrix();
}

/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */
void display() {
    // glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on Model-View matrix
    glLoadIdentity();                       // Reset the model-view matrix

    // default arguments of gluLookAt
    // gluLookAt(0,0,0, 0,0,-100, 0,1,0);

    // control viewing (or camera)
    gluLookAt(pos.x,pos.y,pos.z,
              pos.x+l.x,pos.y+l.y,pos.z+l.z,
              u.x,u.y,u.z);
    // draw
    glScaled(4,4,4);
    glRotated(rotateAngle,0,1,0);
    if(isAxes)
        drawAxes();
    drawSphere();
    drawOctahedron();
    drawCylinder();
    
    glutSwapBuffers();  // Render now
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshapeListener(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}

/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y) {
    double rate = 0.01;
    switch (key) {
        //Octahedron to Sphere
        case ',': 
            octScale = max(0.0, octScale-1.0/STEP_COUNT);
            octTranslate = min(1.0/3.0, octTranslate+1.0/(3.0*STEP_COUNT));
            sphereScale = min(1.0 , sphereScale+1.0/STEP_COUNT);
            sphereTranslate = max(0.0, sphereTranslate-1.0/STEP_COUNT);
            cylinderTranslate = max(0.0, cylinderTranslate-0.5/STEP_COUNT);
            break;
        //Sphere to Octahedron
        case '.': 
            octScale = min(1.0, octScale+1.0/STEP_COUNT);
            octTranslate = max(0.0, octTranslate-1.0/(3.0*STEP_COUNT));
            sphereScale = max(0.0 , sphereScale-1.0/STEP_COUNT);
            sphereTranslate = min(1.0, sphereTranslate+1.0/STEP_COUNT);
            cylinderTranslate = min(0.5, cylinderTranslate+0.5/STEP_COUNT);
            break;

        //control translation
        case '1':
            r.x = r.x*cos(rate)+l.x*sin(rate);
            r.y = r.y*cos(rate)+l.y*sin(rate);
            r.z = r.z*cos(rate)+l.z*sin(rate);

            l.x = l.x*cos(rate)-r.x*sin(rate);
            l.y = l.y*cos(rate)-r.y*sin(rate);
            l.z = l.z*cos(rate)-r.z*sin(rate);
            break;

        case '2':
            r.x = r.x*cos(-rate)+l.x*sin(-rate);
            r.y = r.y*cos(-rate)+l.y*sin(-rate);
            r.z = r.z*cos(-rate)+l.z*sin(-rate);

            l.x = l.x*cos(-rate)-r.x*sin(-rate);
            l.y = l.y*cos(-rate)-r.y*sin(-rate);
            l.z = l.z*cos(-rate)-r.z*sin(-rate);
            break;

        case '3':
            l.x = l.x*cos(rate)+u.x*sin(rate);
            l.y = l.y*cos(rate)+u.y*sin(rate);
            l.z = l.z*cos(rate)+u.z*sin(rate);

            u.x = u.x*cos(rate)-l.x*sin(rate);
            u.y = u.y*cos(rate)-l.y*sin(rate);
            u.z = u.z*cos(rate)-l.z*sin(rate);
            break;

        case '4':
            l.x = l.x*cos(-rate)+u.x*sin(-rate);
            l.y = l.y*cos(-rate)+u.y*sin(-rate);
            l.z = l.z*cos(-rate)+u.z*sin(-rate);

            u.x = u.x*cos(-rate)-l.x*sin(-rate);
            u.y = u.y*cos(-rate)-l.y*sin(-rate);
            u.z = u.z*cos(-rate)-l.z*sin(-rate);
            break;

        case '5':
            u.x = u.x*cos(rate)+r.x*sin(rate);
            u.y = u.y*cos(rate)+r.y*sin(rate);
            u.z = u.z*cos(rate)+r.z*sin(rate);

            r.x = r.x*cos(rate)-u.x*sin(rate);
            r.y = r.y*cos(rate)-u.y*sin(rate);
            r.z = r.z*cos(rate)-u.z*sin(rate);
            break;

        case '6':
            u.x = u.x*cos(-rate)+r.x*sin(-rate);
            u.y = u.y*cos(-rate)+r.y*sin(-rate);
            u.z = u.z*cos(-rate)+r.z*sin(-rate);

            r.x = r.x*cos(-rate)-u.x*sin(-rate);
            r.y = r.y*cos(-rate)-u.y*sin(-rate);
            r.z = r.z*cos(-rate)-u.z*sin(-rate);
            break;

        
        case 'a':
            rotateAngle-=10;   
            break;

        case 'd':
            rotateAngle+=10;   
            break;

        // Control what is shown
        case 'A':
            isAxes = !isAxes;   // show/hide Axes if 'a' is pressed
            break;

        // Control exit
        case 27:    // ESC key
            exit(0);    // Exit window
            break;
    }
    glutPostRedisplay();    // Post a paint request to activate display()
}

/* Callback handler for special-key event */
void specialKeyListener(int key, int x,int y) {
    switch (key) {
        //control rotation
        case GLUT_KEY_UP:		//down arrow key
			pos.x+=l.x;
			pos.y+=l.y;
			pos.z+=l.z;
			break;
		case GLUT_KEY_DOWN:		// up arrow key
			pos.x-=l.x;
			pos.y-=l.y;
			pos.z-=l.z;
			break;

		case GLUT_KEY_RIGHT:
			pos.x+=r.x;
			pos.y+=r.y;
			pos.z+=r.z;
			break;
		case GLUT_KEY_LEFT :
			pos.x-=r.x;
			pos.y-=r.y;
			pos.z-=r.z;
			break;
		case GLUT_KEY_PAGE_UP:
		    pos.x+=u.x;
			pos.y+=u.y;
			pos.z+=u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
            pos.x-=u.x;
			pos.y-=u.y;
			pos.z-=u.z;
			break;
    }
    glutPostRedisplay();    // Post a paint request to activate display()
}

/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
    pos.x=1;pos.y=1;pos.z=20;

    l.x=0;l.y=0;l.z=-1;
    u.x=0;u.y=1;u.z=0;
    r.x=1;r.y=0;r.z=0;
    
    glutInit(&argc, argv);                      // Initialize GLUT
    glutInitWindowSize(640, 640);               // Set the window's initial width & height
    glutInitWindowPosition(100, 50);            // Position the window's initial top-left corner
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color
    glutCreateWindow("Magic Cube");      // Create a window with the given title
    glutDisplayFunc(display);                   // Register display callback handler for window re-paint
    glutReshapeFunc(reshapeListener);           // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener);         // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);        // Register callback handler for special-key event
    initGL();                                   // Our own OpenGL initialization
    glutMainLoop();                             // Enter the event-processing loop
    return 0;
}

/*
g++ $1.cpp -o demo -lglut -lGLU -lGL
./demo
*/