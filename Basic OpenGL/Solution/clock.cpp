#include<iostream>
#include <ctime>
#include<cmath>
#include<GL/glut.h>
using namespace std;

float hourAngle=0.0, minuteAngle=0.0, secondAngle=0.0, pendulumAngle=0.0;
float hourAngleRate=0.00833, minuteAngleRate=0.1, secondAngleRate=6.0;


float thetaMax = 30.0;
float theta = 0.0;
float g = 9.8, L = 0.99;
long long t = 0;

int timerTime = 25;


void drawLine(float x1, float y1, float x2, float y2){
    glBegin(GL_LINES);
        glVertex2f(x1, y1);
        glVertex2f(x2, y2);
    glEnd();
}

void drawCircle(float cx, float cy , float r){
    glBegin(GL_LINE_LOOP); 
        for (float theta = 0; theta < 360; theta += 10) {
            float x = cx + r * cos(theta/180*M_PI);
            float y = cy + r * sin(theta/180*M_PI);
            glVertex2f(x, y);
        }
    glEnd();
}

void drawFilledCircle(float cx, float cy , float r){
    glBegin(GL_POLYGON);
        for (float theta = 0; theta < 360; theta += 10) {
            float x = cx + r * cos(theta/180*M_PI);
            float y = cy + r * sin(theta/180*M_PI);
            glVertex2f(x, y);
        }
    glEnd();
}

void drawRectangle(float x1, float y1, float x2, float y2){
    glBegin(GL_QUADS);
        glVertex2d(x1, y1);
        glVertex2d(x1, y2);
        glVertex2d(x2, y2);
        glVertex2d(x2, y1);
    glEnd();
}

void drawTriangle(float x1, float y1, float x2, float y2, float x3, float y3 ){
    glBegin(GL_TRIANGLES);
        glVertex2d(x1, y1);
        glVertex2d(x2, y2);
        glVertex2d(x3, y3);
    glEnd();
}

void drawDesign(){
    //Bird house
    glBegin(GL_POLYGON);
        glColor3f(0.5,0.0,0.0);
        glVertex2d(-0.4, -0.4);
        glVertex2d(-0.5, 0.2);
        glVertex2d(0, 0.5);
        glColor3f(0.4, 0.0, 0.0);
        glVertex2d(0.5, 0.2);
        glVertex2d(0.4, -0.4);
    glEnd();

    //Left hood of the house
    glBegin(GL_POLYGON);
        glColor3f(0.2, 0.1, 0.0);
        glVertex2d(0, 0.5);
        glVertex2d(0, 0.55);
        glColor3f(0.25, 0.15, 0.05);
        glVertex2d(-0.6, 0.2);
        glVertex2d(-0.5, 0.2);
    glEnd();

    //Right hood of the house
    glBegin(GL_POLYGON);
        glColor3f(0.2, 0.1, 0.0);
        glVertex2d(0, 0.5);
        glVertex2d(0, 0.55);
        glColor3f(0.15, 0.05, 0.0);
        glVertex2d(0.6, 0.2);
        glVertex2d(0.5, 0.2);
    glEnd();

    //Border of the bird house
    glLineWidth(0.5);
    glBegin(GL_LINE_LOOP);
        glColor3f(0.15, 0.05, 0.0);
        glVertex2d(-0.4, -0.4);
        glVertex2d(-0.5, 0.2);
        glVertex2d(0, 0.5);
        glVertex2d(0.5, 0.2);
        glVertex2d(0.4, -0.4);
    glEnd();

    //Floor
    glLineWidth(0.5);
    glBegin(GL_POLYGON);
        glColor3f(0.25, 0.15, 0.05);
        glVertex2d(-0.42, -0.4);
        glVertex2d(-0.45, -0.5);
        glColor3f(0.15, 0.05, 0.0);
        glVertex2d(0.42, -0.4);
        glVertex2d(0.45, -0.5);
    glEnd();
    
}

void display(){
    glClearColor(1, 0.97, 0.97, 1.0); 
    glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glMatrixMode(GL_MODELVIEW);             // To operate on Model-View matrix
    glLoadIdentity();                       // Reset the model-view matrix

    glTranslated(0, 0.3, 0);
    drawDesign();
    

    //clock
    glColor3f(1, 1, 1);
    drawFilledCircle(0, 0, 0.3);

    glColor3f(0.15, 0.05, 0.0);
    glLineWidth(3);

    drawCircle(0, 0, 0.3);
    
    //clock hands
    for(int i=0;i<360;i+=30){
        glPushMatrix();
        glRotated(i, 0,0, 1);
        glTranslated(0, 0.3, 0);
        if(i%90==0) drawLine(0, 0, 0,-0.05);
        else drawLine(0, 0, 0,-0.03);
        glPopMatrix();
    }

    //clock center
    drawRectangle(-0.015, -0.015, 0.015, 0.015);

    //hour
    glPushMatrix();
    
    glRotated(hourAngle, 0 ,0 ,1);
    drawRectangle(-0.01, 0, 0.01, 0.1);
    drawTriangle(-0.01, 0.1, 0.01, 0.1,0,0.17);

    glPopMatrix();
    
    //minute
    glPushMatrix();

    glRotated(minuteAngle, 0 ,0 ,1);
    drawRectangle(-0.007, 0, 0.007, 0.15);
    drawTriangle(-0.007, 0.15, 0.007, 0.15,0,0.20);

    glPopMatrix();


    //second
    glPushMatrix();

    glRotated(secondAngle, 0 ,0 ,1);
    drawRectangle(-0.003, 0, 0.003, 0.18);
    drawTriangle(-0.003, 0.18, 0.003, 0.18,0,0.23);

    glPopMatrix();

    //pendulum
    glPushMatrix();

    glTranslated(0, -0.42, 0);
    glRotated(theta, 0, 0, 1);

    drawFilledCircle(0, 0, 0.03);
    drawRectangle(-0.02, 0, 0.02, -0.45);
    drawFilledCircle(0, -0.45, 0.1);

    glPopMatrix();


    glFlush();

}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }
}

/* Called back when timer expired */
void timer(int value) {
    glutTimerFunc(timerTime, timer, 0); // Call next 'timer' milliseconds later
    glutPostRedisplay();    // Post re-paint request to activate display()
    if(t%1000==0){
        hourAngle -= hourAngleRate;
        minuteAngle -= minuteAngleRate;
        secondAngle -= secondAngleRate;
    }
    
    float omega = sqrt(g/L); // T = 2*pi*sqrt(L/g), omega = sqrt(g/L)
    theta = thetaMax * cos(omega * t/1000.0);
    t = t+timerTime;
}


int main(int argc, char** argv){

    //get current time
    time_t now = time(0);
    char* date_time = ctime(&now);

    int hour = (date_time[11]-'0')*10+(date_time[12]-'0');
    if(hour >= 12) hour -= 12;
    int minute = (date_time[14]-'0')*10+(date_time[15]-'0');
    int second = (date_time[17]-'0')*10+(date_time[18]-'0');

    secondAngle -= 6.0*second;
    minuteAngle -= 0.1*(second + minute*60);
    hourAngle -= 0.00833 * (second + minute*60 + hour*3600);

    //window config
    glutInit(&argc, argv);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(100,50);
    glutCreateWindow("Clock");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(0, timer, 0); 
    glutMainLoop();
}