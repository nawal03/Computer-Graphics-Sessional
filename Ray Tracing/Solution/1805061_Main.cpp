#include <bits/stdc++.h>
#include <GL/glut.h> // GLUT, include glu.h and gl.h
#include "1805061_Header.h"
#include "1805061_bitmap_image.hpp"
using namespace std;

// Global variables
Vector pos; // position of the eye
Vector l;   // look/forward direction
Vector r;   // right direction
Vector u;   // up direction

int near_distance, far_distance;
double fovz;
double aspectratio;
int recur_level;
int npix;
int screen_height, screen_width;
bool isAxes = true;
extern vector<Object *> objects;
extern vector<PointLight *> pointlights;
Floor *flr;

void capture()
{
    bitmap_image *image;
    image = new bitmap_image(npix, npix);
    for (int i = 0; i < npix; i++)
    {
        for (int j = 0; j < npix; j++)
            image->set_pixel(i, j, 0, 0, 0);
    }

    // window screen and image pixel mapping
    // npix -> number of pixels in the image
    double near_distance = (screen_height / 2.0) / tan(fovz * M_PI / 360.0);

    double dx = (1.0 * screen_width / npix) * aspectratio;
    double dz = (1.0 * screen_height / npix);

    Vector mid_point = pos + l * near_distance;
    Vector top_left = mid_point - r * dx * (npix / 2) + u * dz * (npix / 2);
    top_left = top_left + r * (dx / 2.0) - u * (dz / 2.0); // middle point of a pixel

    for (int x = 0; x < npix; x++)
    {
        for (int z = 0; z < npix; z++)
        {
            Vector point = top_left + r * dx * x - u * dz * z;
            Ray ray(pos, point - pos);

            double t_min = inf;
            Color color_min(0, 0, 0);

            for (auto obj : objects)
            {
                Color color(0, 0, 0);
                double t = obj->intersect(ray, color, recur_level);
                if (t >= 0 && t < t_min)
                {
                    t_min = t;
                    color_min = color;
                }
            }

            image->set_pixel(x, z, color_min.r * 255, color_min.g * 255, color_min.b * 255);
        }
    }

    image->save_image("out.bmp");
    delete image;
}

void loadData()
{
    ifstream in("description.txt");

    in >> near_distance >> far_distance >> fovz >> aspectratio >> recur_level >> npix;

    double width, ambient, diffuse, specular, reflection, shine;

    in >> width >> ambient >> diffuse >> reflection;
    flr = new Floor(500 * width, width, ambient, diffuse, reflection);
    objects.push_back(flr);

    int nobj;
    in >> nobj;
    while (nobj--)
    {
        string tobj;
        in >> tobj;
        if (tobj == "sphere")
        {
            Vector center;
            double radius;
            Color color;
            in >> center.x >> center.y >> center.z;
            in >> radius;
            in >> color.r >> color.g >> color.b;
            in >> ambient >> diffuse >> specular >> reflection >> shine;
            objects.push_back(new Sphere(center, radius, color, ambient, diffuse, specular, reflection, shine));
        }
        else if (tobj == "cube")
        {
            Vector bottom_left;
            double side;
            Color color;
            in >> bottom_left.x >> bottom_left.y >> bottom_left.z;
            in >> side;
            in >> color.r >> color.g >> color.b;
            in >> ambient >> diffuse >> specular >> reflection >> shine;
            objects.push_back(new Cube(bottom_left, side, color, ambient, diffuse, specular, reflection, shine));
        }
        else if (tobj == "pyramid")
        {
            Vector lowest_point;
            double width, height;
            Color color;
            in >> lowest_point.x >> lowest_point.y >> lowest_point.z;
            in >> width >> height;
            in >> color.r >> color.g >> color.b;
            in >> ambient >> diffuse >> specular >> reflection >> shine;
            objects.push_back(new Pyramid(lowest_point, width, height, color, ambient, diffuse, specular, reflection, shine));
        }
    }

    int npl;
    in >> npl;
    while (npl--)
    {
        Vector position;
        double falloff;
        in >> position.x >> position.y >> position.z >> falloff;
        pointlights.push_back(new PointLight(position, falloff));
    }

    int nsl;
    in >> nsl;
    while (nsl--)
    {
        Vector position;
        double falloff;
        Vector light_direction;
        double cutoff_angle;
        in >> position.x >> position.y >> position.z >> falloff;
        in >> light_direction.x >> light_direction.y >> light_direction.z;
        in >> cutoff_angle;
        pointlights.push_back(new SpotLight(position, falloff, light_direction, cutoff_angle));
    }

    in.close();
}

/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes()
{
    glLineWidth(1);
    glBegin(GL_LINES);
    glColor3f(1, 0, 0); // Red
    // X axis
    glVertex3f(1000, 0, 0);
    glVertex3f(-1000, 0, 0);

    glColor3f(0, 1, 0); // Green
    // Y axis
    glVertex3f(0, 1000, 0);
    glVertex3f(0, -1000, 0);

    glColor3f(0, 0, 1); // Blue
    // Z axis
    glVertex3f(0, 0, 1000);
    glVertex3f(0, 0, -1000);
    glEnd();
}

void drawObjects()
{
    for (auto obj : objects)
    {
        obj->draw();
    }

    for (auto pl : pointlights)
    {
        pl->draw();
    }
}

/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW); // To operate on Model-View matrix
    glLoadIdentity();           // Reset the model-view matrix

    // control viewing (or camera)
    gluLookAt(pos.x, pos.y, pos.z,
              pos.x + l.x, pos.y + l.y, pos.z + l.z,
              u.x, u.y, u.z);

    // draw
    if (isAxes)
        drawAxes();
    drawObjects();

    glutSwapBuffers(); // Render now
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshapeListener(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0)
        height = 1; // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity();            // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(fovz, aspectratio, near_distance, far_distance);
}

/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y)
{
    double rate = 0.01;
    switch (key)
    {
    case ' ':
        flr->isTexture = !flr->isTexture;
        break;
    case '0':
        cout << "Capturing..." << endl;
        capture();
        cout << "Done" << endl;
        break;
    case '1':
        // u constant
        l = l * cos(rate) + u.cross(l) * sin(rate);
        r = r * cos(rate) + u.cross(r) * sin(rate);
        break;
    case '2':
        // u constant
        l = l * cos(-rate) + u.cross(l) * sin(-rate);
        r = r * cos(-rate) + u.cross(r) * sin(-rate);
        break;
    case '3':
        // r constant
        u = u * cos(rate) + r.cross(u) * sin(rate);
        l = l * cos(rate) + r.cross(l) * sin(rate);
        break;
    case '4':
        // r constant
        u = u * cos(-rate) + r.cross(u) * sin(-rate);
        l = l * cos(-rate) + r.cross(l) * sin(-rate);
        break;
    case '5':
        // l constant
        r = r * cos(rate) + l.cross(r) * sin(rate);
        u = u * cos(rate) + l.cross(u) * sin(rate);
        break;
    case '6':
        // l constant
        r = r * cos(-rate) + l.cross(r) * sin(-rate);
        u = u * cos(-rate) + l.cross(u) * sin(-rate);
        break;
    case 'A':
        isAxes = !isAxes; // show/hide Axes if 'a' is pressed
        break;
    default:
        break;
    }
    glutPostRedisplay(); // Post a paint request to activate display()
}

/* Callback handler for special-key event */
void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_DOWN:
        pos = pos - l;
        break;
    case GLUT_KEY_UP:
        pos = pos + l;
        break;
    case GLUT_KEY_RIGHT:
        pos = pos + r;
        break;
    case GLUT_KEY_LEFT:
        pos = pos - r;
        break;
    case GLUT_KEY_PAGE_UP:
        pos = pos + u;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos - u;
        break;
    default:
        break;
    }
    glutPostRedisplay(); // Post a paint request to activate display()
}

/* Callback handler for mouse click event */
void mouseListener(int button, int state, int x, int y)
{
    switch (button)
    {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN)
        {
            isAxes = !isAxes;
        }
        break;
    default:
        break;
    }
    glutPostRedisplay(); // Post a paint request to activate display()
}

void initGL()
{
    pos = Vector(0.0, -150.0, 50.0);
    r = Vector(1, 0, 0);
    u = Vector(0, 0, 1);
    l = Vector(0, 1, 0);
    loadData();
}

void exitGL()
{
    for (auto obj : objects)
        delete obj;
    for (auto pl : pointlights)
        delete pl;
    objects.clear();
    pointlights.clear();
}

int main(int argc, char **argv)
{
    initGL();
    glutInit(&argc, argv);

    screen_height = npix, screen_width = npix;
    glutInitWindowSize(screen_height, screen_width);
    glutInitWindowPosition(100, 50);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color
    glutCreateWindow("Ray Tracing");

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling

    glutDisplayFunc(display);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);
    glutReshapeFunc(reshapeListener);
    glutMainLoop(); // The main loop of OpenGL
    atexit(exitGL);
    return 0;
}

/*
g++ $1.cpp -o demo -lglut -lGLU -lGL
./demo
*/
