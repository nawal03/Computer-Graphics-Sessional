#include <bits/stdc++.h>
#include <GL/glut.h>
#define STB_IMAGE_IMPLEMENTATION
#include "1805061_stb_image.h"
#define inf 1e8
#define epsi 1e-8
using namespace std;

class Vector{
public:
    double x, y, z;
    Vector(double x=0, double y=0, double z=0){
        this->x =x;
        this->y =y;
        this->z =z;
    }

    Vector operator+(Vector const& obj)
    {
        Vector res;
        res.x = x+obj.x;
        res.y = y+obj.y;
        res.z = z+obj.z;
        return res;
    }

    Vector operator-(Vector const& obj)
    {
        Vector res;
        res.x = x-obj.x;
        res.y = y-obj.y;
        res.z = z-obj.z;
        return res;
    }

    Vector operator-()
    {
        Vector res;
        res.x = -x;
        res.y = -y;
        res.z = -z;
        return res;
    }

    Vector operator*(double d)
    {
        Vector res(x*d, y*d, z*d);
        return res;
    }

    double dot(Vector const& obj)
    {
        return x*obj.x + y*obj.y + z*obj.z;
    }

    Vector cross(Vector const& obj){
        Vector res;
        res.x = y*obj.z - z*obj.y;
        res.y = z*obj.x - x*obj.z;
        res.z = x*obj.y - y*obj.x;
        return res;
    }

    double det()
    {
        return sqrt(x*x+y*y+z*z);
    }

    void normalize(){
        double l = det();
        x = x/l;
        y = y/l;
        z = z/l;
    }

    double distance(Vector v){
        return sqrt((v.x-x)*(v.x-x)+(v.y-y)*(v.y-y)+(v.z-z)*(v.z-z));
    }

};

class Color
{
public:
    double r, g, b;

    Color(double r = 0, double g = 0, double b = 0)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    Color operator+(Color const& obj)
    {
        Color res;
        res.r= r+obj.r;
        res.g = g+obj.g;
        res.b= b+obj.b;
        return res;
    }

    Color operator*(double d)
    {
        Color res(r*d, g*d, b*d);
        return res;
    }

    void clamp()
    {
        this->r = min(this->r, 1.0);
        this->r = max(this->r, 0.0);
        this->g = min(this->g, 1.0);
        this->g = max(this->g, 0.0);
        this->b = min(this->b, 1.0);
        this->b = max(this->b, 0.0);
    }

};

class Ray
{
public:
    Vector start;
    Vector dir;

    Ray(Vector start, Vector dir)
    {
        this->start = start;
        this->dir = dir;
        this->dir.normalize();
    }
};

class Object
{
public:
    Vector reference_point;
    double height, width, length;
    Color color;
    double ambient, diffuse, specular, reflection;
    double shine;

    Object(){
    }
    
    Object(Vector reference_point, double height, double width, double length, Color color, double ambient, double diffuse, double specular,
           double reflection, double shine)
    {
        this->reference_point = reference_point;
        this->height = height;
        this->width = width;
        this->length = length;
        this->color = color;
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->reflection = reflection;
        this->shine = shine;
    }

    virtual void draw(){};

    virtual void getIntersectionResult(Ray r, double &t, Vector& intersection, Vector& normal, Color& intersection_color){}

    double intersect(Ray r, Color &clr, int level);

    Color phongModel(Ray r, Vector intersection, Vector normal, Color intersecion_color);
};

class Triangle : public Object{
public:
    Vector a, b, c;
    
    Triangle(Vector a, Vector b, Vector c){
        this->a = a;
        this->b = b;
        this->c = c;
    }

    void getIntersectionResult(Ray r, double &t, Vector& intersection, Vector& normal, Color& intersection_color){
        double d = (a.x - b.x) * ((a.y - c.y) * r.dir.z - (a.z-c.z) * r.dir.y)
                + (a.y - b.y) * ((a.z - c.z) * r.dir.x - (a.x - c.x) * r.dir.z)
                + (a.z - b.z) * ((a.x - c.x) * r.dir.y - (a.y - c.y) * r.dir.x);

        double d_beta = (a.x - r.start.x) * ((a.y - c.y) * r.dir.z - (a.z-c.z) * r.dir.y)
                    + (a.y - r.start.y) * ( (a.z - c.z) * r.dir.x - (a.x - c.x) * r.dir.z)
                    + (a.z - r.start.z) * ((a.x - c.x) * r.dir.y - (a.y - c.y) * r.dir.x);

        double d_gamma = (a.x - b.x) * ((a.y - r.start.y) * r.dir.z - (a.z-r.start.z) * r.dir.y)
                    + (a.y - b.y) * ((a.z - r.start.z) * r.dir.x - (a.x - r.start.x) * r.dir.z)
                    + (a.z - b.z) * ((a.x - r.start.x) * r.dir.y - (a.y - r.start.y) * r.dir.x);

        double d_t = (a.x - b.x) * ((a.y - c.y) * (a.z - r.start.z) - (a.z-c.z) * (a.y - r.start.y))
                + (a.y - b.y) * ((a.z - c.z) * (a.x - r.start.x) - (a.x - c.x) * (a.z - r.start.z))
                + (a.z - b.z) * ((a.x - c.x) * (a.y - r.start.y) - (a.y - c.y) * (a.x - r.start.x));

        if(abs(d)<epsi){
            t = inf+5;
            return;
        }

        double beta = d_beta/d;
        double gamma = d_gamma/d;
        t = d_t/d;

        if(beta+gamma<1 && beta > 0.0 && gamma > 0.0 && t > 0.0)
        {
            intersection = r.start + r.dir * t;
            normal = (b-a).cross(c-a);
            normal.normalize();
            return;
        }

        t = inf + 5;
        return;
    }

};

class Sphere : public Object
{
public:
    Sphere(Vector center, double radius, Color color, double ambient, double diffuse, double specular,
           double reflection, double shine) : Object(center, radius, radius, radius, color, ambient, diffuse, specular, reflection, shine)
    {
    }

    void draw()
    {
        int stacks = 100, slices = 100;
        Vector points[stacks + 1][slices + 1];
        for (int j = 0; j <= stacks; j++)
        {
            double phi = -M_PI / 2.0 + j * M_PI / stacks;
            double r = length * cos(phi);
            double h = length * sin(phi);
            for (int i = 0; i < slices + 1; i++)
            {
                double theta = i * 2.0 * M_PI / slices;
                points[j][i].x = r * cos(theta);
                points[j][i].y = r * sin(theta);
                points[j][i].z = h;
            }
        }

        glColor3f(color.r, color.g, color.b);
        glPushMatrix();
        glTranslated(reference_point.x, reference_point.y, reference_point.z);
        glBegin(GL_QUADS);
        for (int j = 0; j < stacks; j++)
        {
            for (int i = 0; i < slices; i++)
            {
                glVertex3f(points[j][i].x, points[j][i].y, points[j][i].z);
                glVertex3f(points[j][i + 1].x, points[j][i + 1].y, points[j][i + 1].z);
                glVertex3f(points[j + 1][i + 1].x, points[j + 1][i + 1].y, points[j + 1][i + 1].z);
                glVertex3f(points[j + 1][i].x, points[j + 1][i].y, points[j + 1][i].z);
            }
        }
        glEnd();
        glPopMatrix();
    }

    void getIntersectionResult(Ray r, double &t, Vector& intersection, Vector& normal, Color& intersection_color){
        Vector r0 = r.start - reference_point;
        double tp = -r0.dot(r.dir);
        if(tp<0){
            t = inf+5;
            return;
        }

        double r0_r0 = r0.dot(r0);
        double d_d = r0_r0 - tp*tp;
        if(d_d > length * length){
            t = inf+5;
            return;
        }

        double t_prime = sqrt(length * length - d_d);
        if(r0_r0 < length * length) t = tp + t_prime;
        else t = tp - t_prime;
         
        intersection = r.start + r.dir * t;
        normal = intersection - reference_point;
        normal.normalize();
        intersection_color = color;
    }
};

class Pyramid : public Object
{
public:
    Pyramid(Vector lowest_point, double width, double height, Color color, double ambient, double diffuse, double specular,
            double reflection, double shine) : Object(lowest_point, height, width, width, color, ambient, diffuse, specular, reflection, shine)
    {
    }

    void draw()
    {
        glColor3f(color.r, color.g, color.b);
        glPushMatrix();
        glTranslated(reference_point.x, reference_point.y, reference_point.z);
        glBegin(GL_TRIANGLES);
        glVertex3f(0, 0, 0);
        glVertex3f(0, width, 0);
        glVertex3f(width/2, width/2, height);

        glVertex3f(0, 0, 0);
        glVertex3f(width, 0, 0);
        glVertex3f(width/2, width/2, height);

        glVertex3f(width, width, 0);
        glVertex3f(0, width, 0);
        glVertex3f(width/2, width/2, height);

        glVertex3f(width, width, 0);
        glVertex3f(width, 0, 0);
        glVertex3f(width/2, width/2, height);
        glEnd();

        glBegin(GL_QUADS);
        glVertex3f(0, 0, 0);
        glVertex3f(0, width, 0);
        glVertex3f(width, width, 0);
        glVertex3f(width, 0, 0);
        glEnd();
        glPopMatrix();
    }

    void getIntersectionResult(Ray r, double &t, Vector& intersection, Vector& normal, Color& intersection_color){
        t = inf+5;
        double tr;
        Vector intersectionr;
        Vector normalr;
        Color clr;
        
        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(0, width, 0), 
            reference_point+Vector(width/2, width/2, height)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point+Vector(width, 0, 0), reference_point+Vector(width, width, 0), 
            reference_point+Vector(width/2, width/2, height)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(width, 0, 0), 
            reference_point+Vector(width/2, width/2, height)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point+Vector(0, width, 0), reference_point+Vector(width, width, 0), 
            reference_point+Vector(width/2, width/2, height)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(0, width, 0), 
            reference_point+Vector(width, width, 0)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(width, 0, 0), 
            reference_point+Vector(width, width, 0)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }
        intersection_color = color;
    }
};

class Cube : public Object
{
public:
    Cube(Vector bottom_left, double side, Color color, double ambient, double diffuse, double specular,
         double reflection, double shine) : Object(bottom_left, side, side, side, color, ambient, diffuse, specular, reflection, shine)
    {
    }

    void draw()
    {
        glColor3f(color.r, color.g, color.b);
        glPushMatrix();
        glTranslated(reference_point.x, reference_point.y, reference_point.z);
        glBegin(GL_QUADS);
        glVertex3f(0, 0, 0);
        glVertex3f(0, length, 0);
        glVertex3f(length, length, 0);
        glVertex3f(length, 0, 0);

        glVertex3f(0, 0, 0);
        glVertex3f(0, length, 0);
        glVertex3f(0, length, length);
        glVertex3f(0, 0, length);

        glVertex3f(0, 0, 0);
        glVertex3f(0, 0, length);
        glVertex3f(length, 0, length);
        glVertex3f(length, 0, 0);
        glEnd();

        glTranslated(length, length, length);
        glBegin(GL_QUADS);
        glVertex3f(0, 0, 0);
        glVertex3f(0, -length, 0);
        glVertex3f(-length, -length, 0);
        glVertex3f(-length, 0, 0);

        glVertex3f(0, 0, 0);
        glVertex3f(0, -length, 0);
        glVertex3f(0, -length, -length);
        glVertex3f(0, 0, -length);

        glVertex3f(0, 0, 0);
        glVertex3f(0, 0, -length);
        glVertex3f(-length, 0, -length);
        glVertex3f(-length, 0, 0);
        glEnd();
        glPopMatrix();
    }

    void getIntersectionResult(Ray r, double &t, Vector& intersection, Vector& normal, Color& intersection_color){
        t = inf+5;
        double tr;
        Vector intersectionr;
        Vector normalr;
        Color clr;
        
        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(length, 0, 0), 
            reference_point+Vector(length, length, 0)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(0, length, 0), 
            reference_point+Vector(length, length, 0)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }


        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(0, length, 0), 
            reference_point+Vector(0, length, length)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(0, 0, length), 
            reference_point+Vector(0, length, length)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }


        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(0, 0, length), 
            reference_point+Vector(length, 0, length)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point+Vector(0, 0, 0), reference_point+Vector(length, 0, 0), 
            reference_point+Vector(length, 0, length)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        
        Vector reference_point2 = reference_point + Vector(length, length, length);

        Triangle(reference_point2+Vector(0, 0, 0), reference_point2+Vector(-length, 0, 0), 
            reference_point2+Vector(-length, -length, 0)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point2+Vector(0, 0, 0), reference_point2+Vector(0, -length, 0), 
            reference_point2+Vector(-length, -length, 0)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }


        Triangle(reference_point2+Vector(0, 0, 0), reference_point2+Vector(0, -length, 0), 
            reference_point2+Vector(0, -length, -length)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point2+Vector(0, 0, 0), reference_point2+Vector(0, 0, -length), 
            reference_point2+Vector(0, -length, -length)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }


        Triangle(reference_point2+Vector(0, 0, 0), reference_point2+Vector(0, 0, -length), 
            reference_point2+Vector(-length, 0, -length)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        Triangle(reference_point2+Vector(0, 0, 0), reference_point2+Vector(-length, 0, 0), 
            reference_point2+Vector(-length, 0, -length)).getIntersectionResult(r, tr, intersectionr, normalr, clr);
        if(tr<t){
            t = tr;
            intersection = intersectionr;
            normal = normalr;
        }

        intersection_color = color;
    }
};

class Floor : public Object
{
public:
    vector< vector<Color> > colorw, colorb, textw, textb;
    bool isTexture;
    double floorWidth, tileWidth;

    Floor(double floorWidth, double tileWidth, double ambient, double diffuse, double reflection) : Object(Vector(-floorWidth / 2, -floorWidth / 2, 0), 0, 0, 0, Color(0, 0, 0), ambient, diffuse, 0, reflection, 0)
    {
        this->isTexture = false;
        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;

        //loading color buffers
        colorw.resize(1);
        colorw[0].push_back(Color(1.0, 1.0, 1.0));
        colorb.resize(1);
        colorb[0].push_back(Color(0.0, 0.0, 0.0));

        int width, height, channels;
        unsigned char *imageData = stbi_load("texture_w.bmp", &width, &height, &channels, STBI_rgb);
        
        // Create a 2D array to hold the RGB pixel values
        textw.resize(height);
        for(int i=0;i<height;i++)
            textw[i].resize(width);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = (i * width + j) * channels;
                textw[i][j].r = imageData[index]/255.0;       // Red channel
                textw[i][j].g = imageData[index + 1]/255.0;   // Green channel
                textw[i][j].b  = imageData[index + 2]/255.0;   // Blue channel
            }
        }
        stbi_image_free(imageData);
        imageData = stbi_load("texture_b.bmp", &width, &height, &channels, STBI_rgb);
    
        // Create a 2D array to hold the RGB pixel values
        textb.resize(height);
        for(int i=0;i<height;i++)
            textb[i].resize(width);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = (i * width + j) * channels;
                textb[i][j].r = imageData[index]/255.0;       // Red channel
                textb[i][j].g = imageData[index + 1]/255.0;   // Green channel
                textb[i][j].b  = imageData[index + 2]/255.0;   // Blue channel
            }
        }
        stbi_image_free(imageData);


    }

    void draw()
    {
        if(!isTexture){
            for (int i = -floorWidth / 2, p = 0; i < floorWidth / 2; i += tileWidth, p++)
            {
                for (int j = -floorWidth / 2, q = 0; j < floorWidth / 2; j += tileWidth, q++)
                {
                    if ((p + q) % 2 == 0)
                        glColor3f(colorw[0][0].r, colorw[0][0].g, colorw[0][0].b);
                    else
                        glColor3f(colorb[0][0].r, colorb[0][0].g, colorb[0][0].b);
                    glBegin(GL_QUADS);
                    glVertex3f(i, j, 0);
                    glVertex3f(i, j + tileWidth, 0);
                    glVertex3f(i + tileWidth, j + tileWidth, 0);
                    glVertex3f(i + tileWidth, j, 0);
                    glEnd();
                }
            }
        }
    }

    void getIntersectionResult(Ray r, double &t, Vector& intersection, Vector& normal, Color& intersection_color){
        if(abs(r.dir.z)<epsi) {
            t = inf+5;
            return;
        }
        t = -r.start.z/r.dir.z;
        intersection = r.start + r.dir * t;
        normal = Vector(0, 0, 1);

        //color
        int i = (intersection.x - reference_point.x) / tileWidth;
        if(intersection.x - reference_point.x<0) i--;
        int j = (intersection.y - reference_point.y) / tileWidth;
        if(intersection.y - reference_point.y<0) j--;
        
        if((i+j) % 2 == 0){
            int imageHeight = 1, imageWidth = 1;
            if(isTexture) imageHeight = textw.size(), imageWidth = textw[0].size();
            int pointx = ((intersection.x - reference_point.x) - i * tileWidth)/tileWidth * imageWidth;
            pointx=(pointx+imageWidth)%imageWidth;
            int pointy = ((intersection.y - reference_point.y) - j * tileWidth)/tileWidth * imageHeight;
            pointy=(pointy+imageHeight)%imageHeight;
            isTexture? intersection_color = textw[pointy][pointx] : intersection_color = colorw[0][0];
        }
        else{
            int imageHeight = 1, imageWidth = 1;
            if(isTexture) imageHeight = textb.size(), imageWidth = textb[0].size();
            int pointx = ((intersection.x - reference_point.x) - i * tileWidth)/tileWidth * imageWidth;
            pointx=(pointx+imageWidth)%imageWidth;
            int pointy = ((intersection.y - reference_point.y) - j * tileWidth)/tileWidth * imageHeight;
            pointy=(pointy+imageHeight)%imageHeight;
            isTexture? intersection_color = textb[pointy][pointx] : intersection_color = colorb[0][0];
        }
    }

    ~Floor(){
        for(int i=0;i<colorw.size();i++) colorw[i].clear();
        for(int i=0;i<colorb.size();i++) colorb[i].clear();
        for(int i=0;i<textw.size();i++) textw[i].clear();
        for(int i=0;i<textb.size();i++) textb[i].clear();

        colorw.clear();
        colorb.clear();
        textw.clear();
        textb.clear();
    }
};

class PointLight
{
public:
    Vector light_pos;
    double falloff;

    PointLight(Vector light_pos, double falloff)
    {
        this->light_pos = light_pos;
        this->falloff = falloff;
    }

    virtual void draw()
    {
        glPushMatrix();
        glColor3f(1, 1, 1);
        glPointSize(10);
        glBegin(GL_POINTS);
        glVertex3f(light_pos.x, light_pos.y, light_pos.z);
        glEnd();
        glPopMatrix();
    }

    virtual bool isPathOpen(Object* obj, Vector intersection, Vector normal);
};

class SpotLight : public PointLight
{
public:
    Vector light_direction;
    double cutoff_angle;

    SpotLight(Vector light_pos, double falloff, Vector light_direction, double cutoff_angle) : PointLight(light_pos, falloff)
    {
        this->light_direction = light_direction;
        this->cutoff_angle = cutoff_angle;
    }

    void draw()
    {
        int segments = 100;
        double radius = tan(cutoff_angle * M_PI / 360), height = 1;
        Vector dir_vector = light_direction - light_pos;
        dir_vector.normalize(); // Ensure the direction vector is normalized

        double tempx = radius, tempy = 0;
        double currx, curry;

        // Calculate rotation axis and angle
        Vector rotation_axis = Vector(0, 0, -1).cross(dir_vector);
        double rotation_angle = acos(Vector(0, 0, -1).dot(dir_vector)) * 180.0 / M_PI;

        glPushMatrix();
        glTranslated(light_pos.x, light_pos.y, light_pos.z);
        glRotated(rotation_angle, rotation_axis.x, rotation_axis.y, rotation_axis.z);
        glScaled(15,15,15);
        glBegin(GL_TRIANGLES);
        for (int i = 1; i <= segments; i++)
        {
            double theta = i * 2.0 * M_PI / segments;
            currx = radius * cos(theta);
            curry = radius * sin(theta);
            GLfloat c = (2 + cos(theta)) / 3;
            glColor3f(c, c, c);
            glVertex3f(0, 0, height / 2);
            glVertex3f(currx, curry, -height / 2);
            glVertex3f(tempx, tempy, -height / 2);
            tempx = currx;
            tempy = curry;
        }
        glEnd();
        glPopMatrix();
    } 

    bool isPathOpen(Object* obj, Vector intersection, Vector normal);
};

vector<Object *> objects;
vector<PointLight *> pointlights;

bool PointLight::isPathOpen(Object* obj, Vector intersection, Vector normal){
    Vector to_obj = intersection - light_pos;
    to_obj.normalize();
    Ray L = Ray(light_pos, to_obj);
    if(normal.dot(-L.dir)<0) return false;

    // check if the ray intersects with any object
    Color c;
    double tthis = obj->intersect(L, c, 0);
    for(auto o: objects)
    {
        if(o == obj) continue;
        double tnow = o->intersect(L, c, 0);
        if(tnow > 0.0 && (tthis - tnow)  > epsi) return false;
    }
    
    return true;
}

bool SpotLight::isPathOpen(Object* obj, Vector intersection, Vector normal){
    Vector to_obj = intersection - light_pos;
    to_obj.normalize();
    Ray L = Ray(light_pos, to_obj);
    if(normal.dot(-L.dir)<0)    return false;
    
    // check if the ray intersects with any object
    Color c;
    double tthis = obj->intersect(L, c, 0);
    for(auto o: objects)
    {
        if(o == obj) continue;
        double tnow = o->intersect(L, c, 0);
        if(tnow > 0.0 && (tthis - tnow)  > epsi)   return false;
        
    }
    
    //check if the object in the proper angle
    Vector dir_vector = light_direction - light_pos;
    dir_vector.normalize();
    double angle = acos(to_obj.dot(dir_vector)) * 180.0 / M_PI;
    if(angle > cutoff_angle) return false;

    return true;
}

Color Object::phongModel(Ray r, Vector intersection, Vector normal, Color intersecion_color){
    double lambert = 0, phong = 0;
    for(auto pl : pointlights){
        if(!pl->isPathOpen(this, intersection, normal)) continue;

        Vector to_source = pl->light_pos - intersection;
        to_source.normalize();

        double distance = intersection.distance(pl->light_pos);
        double scaling_factor = pow(M_E,-distance*distance*pl->falloff);
        lambert += to_source.dot(normal)*scaling_factor;

        Vector reflection_dir = normal * normal.dot(-r.dir)*2 + r.dir;
        reflection_dir.normalize();
        phong += pow(reflection_dir.dot(to_source), shine)*scaling_factor;
    }

    Color clr = intersecion_color * ambient;
    clr.clamp();
    clr = clr + intersecion_color * diffuse * lambert;
    clr.clamp();
    clr = clr + intersecion_color * specular * phong;
    clr.clamp();
    return clr;
}

double Object::intersect(Ray r, Color &clr, int level){
    double t;
    Vector intersection;
    Vector normal;
    Color intersection_color;
    getIntersectionResult(r, t, intersection, normal, intersection_color);

    clr = Color(0, 0, 0);
    if(level==0 || t>inf || t<0.0)return t;

    if(normal.dot(-r.dir) < 0.0){
        normal = -normal;
    }
    clr = phongModel(r, intersection, normal, intersection_color);

    Vector reflection_dir = normal * normal.dot(-r.dir)*2 + r.dir;
    reflection_dir.normalize();
    Ray reflected = Ray(intersection + reflection_dir * 0.00001, reflection_dir);
    Color reflected_clr;

    double tmin = inf;
    Object* obj;
    for(auto o: objects){
        double tnow = o->intersect(reflected, reflected_clr, 0);
        if(tnow > 0.0 && (tmin - tnow)  > epsi){
            tmin = tnow;
            obj = o;
        }
    }

    if(tmin<inf){
        obj->intersect(reflected, reflected_clr, level-1);
        clr = clr + reflected_clr * reflection;  
        clr.clamp();
    }

    return t;
}