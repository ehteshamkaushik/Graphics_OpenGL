#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include <fstream>
#include <windows.h>
#include <glut.h>
#include <vector>
#include "bitmap_image.hpp"
using namespace std;
#define pi (2*acos(0.0))

int drawaxes, drawgrid;
double n, f, fovY, ar, fovX;
int lor;
int pixelCount;
double cbWidth, cbAmbient, cbDiffused, cbReflection;
double height, width;
double farPlaneD;

class Point
{
public:
    double x, y, z;
    Point()
    {

    }
    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    void operator = (const Point &p )
    {
        x = p.x;
        y = p.y;
        z = p.z;
    }

    Point operator+(const Point& p)
    {
        Point point;
        point.x = this->x + p.x;
        point.y = this->y + p.y;
        point.z = this->z + p.z;
        return point;
    }
    Point operator-(const Point& p)
    {
        Point point;
        point.x = this->x - p.x;
        point.y = this->y - p.y;
        point.z = this->z - p.z;
        return point;
    }
    Point operator*(const double& m)
    {
        Point point;
        point.x = this->x * m;
        point.y = this->y * m;
        point.z = this->z * m;
        return point;
    }

    Point operator*(const Point& p)
    {
        Point point;

        point.x = (y*p.z) - (z*p.y);
        point.y = (z*p.x) - (x*p.z);
        point.z = (x*p.y) - (y*p.x);
        double v = (sqrt(point.x*point.x+point.y*point.y+point.z*point.z));
        point = point / v;
        return point;
    }

    Point operator/(const double& m)
    {
        Point point;
        point.x = this->x / m;
        point.y = this->y / m;
        point.z = this->z / m;
        return point;
    }
    double dot(Point b)
    {
        return x*b.x + y*b.y + z*b.z;
    }
    void printPoint()
    {
        cout << "Point : " << x << " " << y << " " << z << endl;
    }
};

Point** pointBuffer;
Point midpoint;

void normalize(Point &n)
{
    double v = sqrt(n.x*n.x + n.y*n.y+n.z*n.z);
    n = n / v;
}

void drawSquare(double x, double y, double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f(x, y , 0);
		glVertex3f(x + a,y, 0);
		glVertex3f(x + a,y + a, 0);
		glVertex3f(x, y + a, 0);
	}glEnd();
}

void drawCheckerBoard(double a)
{
    double x, y;
    int color = 1;
    for(x = -100*a; x <= 100*a; x += a)
    {
        for(y = -100*a; y <= 100*a; y += a)
        {
            glColor3f(color, color, color);
            drawSquare(x, y, a);
            color = 1 - color;
        }
    }
}


class Color {
public:
    double r, g, b;
    Color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color() {
    }
};

Color **textureBuffer;
int theight, twidth;
bool textureEnable = false;

void loadTextureImage()
{
    bitmap_image b_img ("texture.bmp");

    theight = b_img.height();
    twidth = b_img.width();
    textureBuffer = new Color* [twidth];
    for (int i = 0; i < twidth; i++) {
        textureBuffer[i] = new Color [theight];
        for (int j = 0; j < theight; j++) {
            unsigned char r, g, b;
            b_img.get_pixel(i, j, r, g, b);
            Color c(r/255.0, g/255.0, b/255.0);
            textureBuffer[i][j] = c;
        }
    }
}

Color backgroud(0, 0, 0);
Color** pixels;

class Coeff
{
public:
    double ambient, diffuse, specular, reflection;
    Coeff(){}
    Coeff(double a, double d, double s, double r)
    {
        this->ambient = a;
        this->diffuse = d;
        this->specular = s;
        this->reflection = r;
    }
};

class Sphere
{
public:
    Point center;
    double radius;
    Color color;
    Coeff coeff;
    double shininess;
    Sphere(){}
    void printSphere()
    {
        cout << "Sphere : \n";
        cout << center.x << " " << center.y << " " << center.z << endl;
        cout << radius << endl;
        cout << color.r << " " << color.g << " " << color.b << endl;
        cout << coeff.ambient << " " << coeff.diffuse << " " << coeff.specular << " " << coeff.reflection << endl;
        cout << shininess << endl;
    }

    void drawSphere(int slices,int stacks)
    {
        Point points[100][100];
        int i,j;
        double h,r;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }glEnd();
            }
        }
    }
};

class Triangle
{
public:
    Point p1, p2, p3;
    Triangle()
    {

    }
    Triangle(Point p1, Point p2, Point p3)
    {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }
};

class Pyramid
{
public:
    Point lowestPoint;
    double width, height;
    Color color;
    Coeff coeff;
    double shininess;
    Pyramid(){}
    vector<Triangle> triangles;
    void printPyramid()
    {
        cout << "Pyramid : \n";
        cout << lowestPoint.x << " " << lowestPoint.y << " " << lowestPoint.z << endl;
        cout << width << " " <<  height << endl;
        cout << color.r << " " << color.g << " " << color.b << endl;
        cout << coeff.ambient << " " << coeff.diffuse << " " << coeff.specular << " " << coeff.reflection << endl;
        cout << shininess << endl;
    }

    void drawPyramid()
    {
        drawSquare(0, 0, width);
        Point top(width/2, width/2, height);
        glBegin(GL_TRIANGLES);{
            glVertex3f(top.x, top.y , top.z);
            glVertex3f(0, 0, 0);
            glVertex3f(width, 0, 0);

            glVertex3f(top.x, top.y , top.z);
            glVertex3f(0, 0, 0);
            glVertex3f(0, width, 0);

            glVertex3f(top.x, top.y , top.z);
            glVertex3f(0, width, 0);
            glVertex3f(width, width, 0);

            glVertex3f(top.x, top.y , top.z);
            glVertex3f(width, width, 0);
            glVertex3f(width, 0, 0);
        }glEnd();
    }
};

class LightSource
{
public:
    Point position;
    double falloff;
    LightSource(){}
    void printLightSource()
    {
        cout << "Light Source :\n";
        cout << position.x << " " << position.y << " " << position.z << endl;
        cout << falloff << endl;
    }

    void drawSphere(double radius, int slices,int stacks)
    {
        Point points[100][100];
        int i,j;
        double h,r;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }glEnd();
            }
        }
    }
};

class SpotLightSource
{
public:
    Point position;
    double falloff;
    Point lookAt;
    double cutoffAngle;
    SpotLightSource(){}
    void printLightSource()
    {
        cout << "Spot Light Source :\n";
        cout << position.x << " " << position.y << " " << position.z << endl;
        cout << falloff << endl;
        cout << lookAt.x << " " << lookAt.y << " " << lookAt.z << endl;
        cout << cutoffAngle << endl;
    }
    void drawSphere(double radius, int slices,int stacks)
    {
        Point points[100][100];
        int i,j;
        double h,r;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }glEnd();
            }
        }
    }
};

double calcDist(Point p, Point q)
{
    double d = sqrt((p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y) +(p.z-q.z)*(p.z-q.z));
    return d;
}

vector<Sphere> spheres;
vector<Pyramid> pyramids;
vector<LightSource> lightSources;
vector<SpotLightSource> spotLightSources;

Point pos, l, r, u;

void generatePointBuffer()
{
    height = 2 * n * tan((fovY/2) * (pi/180));
    width = 2 * n * tan((fovX/2) * (pi/180));
    midpoint = pos + l * n;
    Point lowerLeft = midpoint - r * (int)(pixelCount/2) * (width/pixelCount);
    lowerLeft = lowerLeft - u * (int)(pixelCount/2) * (height/pixelCount);
    for (int i = 0; i < pixelCount; i++)
    {
        for (int j = 0; j < pixelCount; j++)
        {
            pointBuffer[i][j] = lowerLeft + u * i * height/pixelCount + r * j * width/pixelCount;
        }
    }
}

double matrixDet(double matrix[3][3])
{
    return ((matrix[0][0]*(matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1]))-
            (matrix[0][1]*(matrix[1][0]*matrix[2][2] - matrix[2][0]*matrix[1][2]))+
            (matrix[0][2]*(matrix[1][0]*matrix[2][1] - matrix[2][0]*matrix[1][1])));
}

double triangleIntersection(Triangle triangle, Point currentPoint, Point rayVector)
{
    double matA[3][3] = {{triangle.p1.x - triangle.p2.x, triangle.p1.x - triangle.p3.x, rayVector.x},
                            {triangle.p1.y - triangle.p2.y, triangle.p1.y - triangle.p3.y, rayVector.y},
                            {triangle.p1.z - triangle.p2.z, triangle.p1.z - triangle.p3.z, rayVector.z}};
    double matB[3][1] = {{triangle.p1.x - currentPoint.x},
                            {triangle.p1.y - currentPoint.y},
                            {triangle.p1.z - currentPoint.z}};
    double matrixBeta[3][3];
    double matrixGama[3][3];
    double matrixT[3][3];
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            matrixBeta[i][j] = matA[i][j];
            matrixGama[i][j] = matA[i][j];
            matrixT[i][j] = matA[i][j];
        }
    }
    matrixBeta[0][0] = matB[0][0];
    matrixBeta[1][0] = matB[1][0];
    matrixBeta[2][0] = matB[2][0];

    matrixGama[0][1] = matB[0][0];
    matrixGama[1][1] = matB[1][0];
    matrixGama[2][1] = matB[2][0];

    matrixT[0][2] = matB[0][0];
    matrixT[1][2] = matB[1][0];
    matrixT[2][2] = matB[2][0];

    double detA = matrixDet(matA);
    double detBeta = matrixDet(matrixBeta);
    double detGama = matrixDet(matrixGama);
    double detT = matrixDet(matrixT);
    double beta = detBeta/detA;
    double gama = detGama/detA;
    double t = detT/detA;
    if(beta + gama < 1 && beta > 0 && gama > 0)
        return t;
    else return -1;
}

bool checkIntersection(Point currentPoint, Point rayVector, double tmin)
{
    bool hasIntersection = false;
    normalize(rayVector);
    double t = -(currentPoint.z)/(rayVector.z);
    if (t < tmin && t >= 0)
    {
        return true;
    }

    for (int s = 0; s < spheres.size(); s++)
    {
        Sphere currentSphere = spheres[s];
        int flag;
        Point newPoint = currentPoint - currentSphere.center;
        double temp = newPoint.x * newPoint.x + newPoint.y * newPoint.y + newPoint.z * newPoint.z;
        if(temp < currentSphere.radius*currentSphere.radius)
        {
            flag = 0;
        }
        else flag = 1;
        double tp = - (newPoint.x * rayVector.x + newPoint.y * rayVector.y + newPoint.z * rayVector.z);
        if(flag == 1 && tp < 0)
        {
            continue;
        }
        double d_sqrd = temp - tp * tp;
        if(d_sqrd > currentSphere.radius*currentSphere.radius)
            continue;
        double t_temp = currentSphere.radius * currentSphere.radius - d_sqrd;
        if(flag == 1)
        {
            t = tp - sqrt(t_temp);
        }
        else
        {
            t = tp + sqrt(t_temp);
        }
        if(t < tmin)
        {
            return true;
        }
    }

    for (int p = 0; p < pyramids.size(); p++)
    {
        Pyramid currentPyramid = pyramids[p];
        double pyrT = 100000;
        for(int tri = 0; tri < currentPyramid.triangles.size(); tri++)
        {
            Triangle triangle = currentPyramid.triangles[tri];
            double currT = triangleIntersection(triangle, currentPoint, rayVector);
            if(currT >= 0 && currT < pyrT)
            {
                pyrT = currT;
            }
        }
        if(pyrT < tmin)
        {
            return true;
        }
    }

    return false;
}

Color rayTracing(Point currentPoint, Point rayVector, int depth)
{
    Color color = Color(0, 0, 0);
    if(depth == 0)
    {
        //cout << "returning\n";
        return color;
    }
    normalize(rayVector);
    double tmin = 100000;
    double shininess = 0;
    Point normal = Point(0, 0, 0);
    Coeff coeff;
    double t;
    if(rayVector.z != 0)
    {
        t = -(currentPoint.z)/(rayVector.z);
        if (t < tmin && t >= 0)
        {
            coeff.ambient = cbAmbient;
            coeff.diffuse = cbDiffused;
            coeff.specular = 0;
            coeff.reflection = cbReflection;
            tmin = t;
            normal = Point(0, 0, 1);
            Point intersection = currentPoint + rayVector * t;
            if (textureEnable)
            {
                int r = (int)(abs(intersection.x)) % (int)cbWidth;
                int s = (int)(abs(intersection.y)) % (int)cbWidth;
                if (intersection.x < 0)
                    r = cbWidth - r - 1;
                if (intersection.y < 0)
                    s = cbWidth - s - 1;
                int xInt = (int)(theight*r/cbWidth);
                int yInt = (int)(twidth*s/cbWidth);

                color = textureBuffer[xInt][yInt];

            }
            else
            {
                if(intersection.x * intersection.y >= 0)
                {
                    if (((int)(abs(intersection.x/cbWidth)) + (int)(abs(intersection.y/cbWidth))) % 2 == 0)
                    {
                        color = Color(1, 1, 1);
                    }
                }
                else
                {
                    if (((int)(abs(intersection.x/cbWidth)) + (int)(abs(intersection.y/cbWidth))) % 2 == 1)
                    {
                        color = Color(1, 1, 1);
                    }
                }
            }
        }
    }


    for (int s = 0; s < spheres.size(); s++)
    {
        Sphere currentSphere = spheres[s];
        int flag;
        Point newPoint = currentPoint - currentSphere.center;
        double temp = newPoint.x * newPoint.x + newPoint.y * newPoint.y + newPoint.z * newPoint.z;
        if(temp < currentSphere.radius*currentSphere.radius)
        {
            flag = 0;
        }
        else flag = 1;
        double tp = - (newPoint.x * rayVector.x + newPoint.y * rayVector.y + newPoint.z * rayVector.z);
        if(flag == 1 && tp < 0)
        {
            continue;
        }
        double d_sqrd = temp - tp * tp;
        if(d_sqrd > currentSphere.radius*currentSphere.radius)
            continue;
        double t_temp = currentSphere.radius * currentSphere.radius - d_sqrd;
        if(flag == 1)
        {
            t = tp - sqrt(t_temp);
        }
        else
        {
            t = tp + sqrt(t_temp);
        }
        if(t < tmin)
        {
            tmin = t;
            Point si = currentPoint + rayVector * tmin;
            si = si - currentSphere.center;
            normalize(si);
            normal = si;
            color = currentSphere.color;
            shininess = currentSphere.shininess;
            coeff = currentSphere.coeff;
        }
    }

    for (int p = 0; p < pyramids.size(); p++)
    {
        Pyramid currentPyramid = pyramids[p];
        //double pyrT = 100000;
        for(int tri = 0; tri < currentPyramid.triangles.size(); tri++)
        {
            Triangle triangle = currentPyramid.triangles[tri];
            double currT = triangleIntersection(triangle, currentPoint, rayVector);
            if(currT >= 0 && currT < tmin)
            {
                Point pp = (triangle.p2 - triangle.p1) * (triangle.p3 - triangle.p1);
                normalize(pp);
                normal = pp;
                if(tri < 2)
                    normal = Point(0, 0, -1);
                tmin = currT;
                color = currentPyramid.color;
                shininess = currentPyramid.shininess;
                coeff = currentPyramid.coeff;
            }
        }
    }

    double fpi = -(farPlaneD + l.dot(currentPoint))/(l.dot(rayVector));
    if (fpi >= 0 && fpi < tmin)
    {
        return Color(0, 0, 0);
    }

    if(tmin >= 100000)
    {
        return color;
    }


    Point intersectionPoint = currentPoint + rayVector * tmin;
    double lambert, phong;
    lambert = 0;
    phong = 0;
    Point R = rayVector - normal * rayVector.dot(normal) * 2;
    normalize(R);
    for (int l = 0; l < lightSources.size(); l++)
    {
        LightSource currentLightSource = lightSources[l];
        Point toSource = currentLightSource.position - intersectionPoint;
        Point R_temp = intersectionPoint - currentLightSource.position;
        normalize(R_temp);
        normalize(toSource);
        intersectionPoint = intersectionPoint - rayVector * .001;
        double tSource = (currentLightSource.position.x - intersectionPoint.x)/toSource.x;
        bool hasIntersection = checkIntersection(intersectionPoint, toSource, tSource);
        if (hasIntersection)
        {
            continue;
        }
        double dist = calcDist(intersectionPoint, currentLightSource.position);
        double scalingFactor = exp(-dist*dist*currentLightSource.falloff);
        lambert += toSource.dot(normal)*scalingFactor;
        Point R_prime = normal * R_temp.dot(normal) * 2 - R_temp;
        normalize(R_prime);
        phong += pow(R.dot(toSource), shininess)*scalingFactor;
    }

    for (int sl = 0; sl < spotLightSources.size(); sl++)
    {
        SpotLightSource currentSpotLightSource = spotLightSources[sl];
        Point toSource = currentSpotLightSource.position - intersectionPoint;
        Point R_temp = intersectionPoint - currentSpotLightSource.position;
        normalize(R_temp);
        normalize(toSource);
        intersectionPoint = intersectionPoint - rayVector * .001;
        double tSource = (currentSpotLightSource.position.x - intersectionPoint.x)/toSource.x;
        bool hasIntersection = checkIntersection(intersectionPoint, toSource, tmin);
        if (hasIntersection)
        {
            continue;
        }
        Point v1 = intersectionPoint - currentSpotLightSource.position;
        normalize(v1);
        Point v2 = currentSpotLightSource.lookAt;
        normalize(v2);
        double angle = acos(v1.dot(v2));
        angle = angle * (180/pi);
        if (angle > currentSpotLightSource.cutoffAngle)
        {
            continue;
        }
        double dist = calcDist(intersectionPoint, currentSpotLightSource.position);
        double scalingFactor = exp(-dist*dist*currentSpotLightSource.falloff);
        lambert += toSource.dot(normal)*scalingFactor;
        Point R_prime = normal * R_temp.dot(normal) * 2 - R_temp;
        normalize(R_prime);
        phong += pow(R.dot(toSource), shininess)*scalingFactor;
    }
    Point reflectedPoint = intersectionPoint - rayVector * .001;
    //cout << "Calling\n";

    Color reflectedColor = rayTracing(reflectedPoint, R, depth - 1);
    color.r = coeff.ambient*color.r + coeff.diffuse*lambert *color.r + coeff.specular*phong *color.r + coeff.reflection*reflectedColor.r;
    color.g = coeff.ambient*color.g + coeff.diffuse*lambert *color.g + coeff.specular*phong *color.g + coeff.reflection*reflectedColor.g;
    color.b = coeff.ambient*color.b + coeff.diffuse*lambert *color.b + coeff.specular*phong *color.b + coeff.reflection*reflectedColor.b;

    return color;
}

void generateRay()
{
    Point farPoint = pos + l * f;
    farPlaneD = - (farPoint.dot(l));
    for (int i = 0; i < pixelCount; i++)
    {
        if ((i+1) %100 == 0)
        {
            cout << "Rendering : "  << (i+1)*1.0/pixelCount*100 << "%\n";
        }
        for (int j = 0; j < pixelCount; j++)
        {

            Point currentPoint = pointBuffer[i][j];
            Point rayVector = currentPoint - pos;
            Color color = rayTracing(currentPoint, rayVector, 3);
            color.r = color.r*255;
            color.g = color.g*255;
            color.b = color.b*255;
            pixels[j][pixelCount - 1 -i] = color;
        }
    }
}

void generateBitmapImage()
{
    bitmap_image image(pixelCount, pixelCount);
    for (int x = 0; x < pixelCount; x++) {
        for (int y = 0; y < pixelCount; y++) {
            image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
        }
    }
    image.save_image("out.bmp");
}

void rotateVector(Point &l, Point &r, Point &u, double a)
{
    Point temp;
    temp = l*acos(a*pi/180) + u*asin(a*pi/180);
    double v = (sqrt(temp.x*temp.x+temp.y*temp.y+temp.z*temp.z));
    temp = temp / v;
    l = temp;
    u = r * l;
}


void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 500,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        case '1':
		    rotateVector(r, u, l, 3);
			break;
        case '2':
			rotateVector(r, u, l, -3);
			break;
		case '3':
		    rotateVector(l, r, u, 3);
			break;
        case '4':
			rotateVector(l, r, u, -3);
			break;
        case '5':
			rotateVector(u, l, r, 3);
			break;
        case '6':
			rotateVector(u, l, r, -3);
			break;
        case '0':
            generatePointBuffer();
            generateRay();
            generateBitmapImage();
            cout << "Rendering Complete\n";
            break;
        case ' ':
            if(!textureEnable)
            {
                textureEnable = true;
            }
            else
            {
                textureEnable = false;
            }
            break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
            pos = pos - l;
			break;
		case GLUT_KEY_UP:		// up arrow key
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

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	if (!textureEnable)
    {
        drawCheckerBoard(cbWidth);
    }
	for(int i = 0; i < spheres.size(); i++)
    {
        Sphere sphere = spheres[i];
        glPushMatrix();
        {
            glColor3f(sphere.color.r, sphere.color.g, sphere.color.b);
            glTranslatef(sphere.center.x, sphere.center.y, sphere.center.z);
            sphere.drawSphere(20, 20);
        }
        glPopMatrix();
    }

    for(int i = 0; i < pyramids.size(); i++)
    {
        Pyramid pyramid = pyramids[i];
        glPushMatrix();
        {
            glColor3f(pyramid.color.r, pyramid.color.g, pyramid.color.b);
            glTranslatef(pyramid.lowestPoint.x, pyramid.lowestPoint.y, pyramid.lowestPoint.z);
            pyramid.drawPyramid();
        }
        glPopMatrix();
    }

    for(int i = 0; i < lightSources.size(); i++)
    {
        LightSource ls = lightSources[i];
        glPushMatrix();
        {
            glColor3f(1, 1, 1);
            glTranslatef(ls.position.x, ls.position.y, ls.position.z);
            ls.drawSphere(5, 10, 10);
        }
        glPopMatrix();
    }

    for(int i = 0; i < spotLightSources.size(); i++)
    {
        SpotLightSource ls = spotLightSources[i];
        glPushMatrix();
        {
            glColor3f(1, 1, 1);
            glTranslatef(ls.position.x, ls.position.y, ls.position.z);
            ls.drawSphere(5, 10, 10);
        }
        glPopMatrix();
    }

	glutSwapBuffers();
}


void animate(){
	//angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void readInput()
{
    ifstream scene;
    scene.open ("description.txt");
    scene >> n >> f >> fovY >> ar;
    //cout << n << " " << f << " " << fovY << " " << ar << endl;
    fovX = ar*fovY;
    scene >> lor;
    //cout << lor << endl;
    scene >> pixelCount;
    //cout << pixelCount << endl;
    pixels = new Color*[pixelCount];
    pointBuffer = new Point*[pixelCount];
    for (int i = 0; i < pixelCount; i++)
    {
        pixels[i] = new Color [pixelCount];
        for (int j = 0; j < pixelCount; j++)
        {
            pixels[i][j] = backgroud;
        }
        pointBuffer[i] = new Point[pixelCount];
        for (int j = 0; j < pixelCount; j++)
        {
            pointBuffer[i][j] = Point(0, 0, 0);
        }
    }
    scene >> cbWidth;
    //cout << cbWidth << endl;
    scene >> cbAmbient >> cbDiffused >> cbReflection;
    //cout << cbAmbient << " " << cbDiffused << " " << cbReflection << endl;
    double objectsCnt;
    scene >> objectsCnt;
    //cout << objectsCnt << endl;
    for(int i = 0; i < objectsCnt; i++)
    {
        string objName;
        scene >> objName;
        if(objName == "sphere")
        {
            Sphere sphere;
            scene >> sphere.center.x >> sphere.center.y >> sphere.center.z;
            scene >> sphere.radius;
            scene >> sphere.color.r >> sphere.color.g >> sphere.color.b;
            scene >> sphere.coeff.ambient >> sphere.coeff.diffuse >> sphere.coeff.specular >> sphere.coeff.reflection;
            scene >> sphere.shininess;
            //sphere.printSphere();
            spheres.push_back(sphere);
        }
        else if(objName == "pyramid")
        {
            Pyramid pyramid;
            scene >> pyramid.lowestPoint.x >> pyramid.lowestPoint.y >> pyramid.lowestPoint.z;
            scene >> pyramid.width >> pyramid.height;
            scene >> pyramid.color.r >> pyramid.color.g >> pyramid.color.b;
            scene >> pyramid.coeff.ambient >> pyramid.coeff.diffuse >> pyramid.coeff.specular >> pyramid.coeff.reflection;
            scene >> pyramid.shininess;
            Point top = Point(pyramid.lowestPoint.x + pyramid.width/2, pyramid.lowestPoint.y + pyramid.width/2, pyramid.lowestPoint.z + pyramid.height);
            pyramid.triangles.push_back(Triangle(Point(pyramid.lowestPoint.x, pyramid.lowestPoint.y, pyramid.lowestPoint.z),
                                                  Point(pyramid.lowestPoint.x + pyramid.width, pyramid.lowestPoint.y, pyramid.lowestPoint.z),
                                                  Point(pyramid.lowestPoint.x + pyramid.width, pyramid.lowestPoint.y + pyramid.width, pyramid.lowestPoint.z)));
            pyramid.triangles.push_back(Triangle(Point(pyramid.lowestPoint.x + pyramid.width, pyramid.lowestPoint.y + pyramid.width, pyramid.lowestPoint.z),
                                                 Point(pyramid.lowestPoint.x, pyramid.lowestPoint.y + pyramid.width, pyramid.lowestPoint.z),
                                                 Point(pyramid.lowestPoint.x, pyramid.lowestPoint.y, pyramid.lowestPoint.z)));

            pyramid.triangles.push_back(Triangle(top, Point(pyramid.lowestPoint.x, pyramid.lowestPoint.y, pyramid.lowestPoint.z),
                                                  Point(pyramid.lowestPoint.x + pyramid.width, pyramid.lowestPoint.y, pyramid.lowestPoint.z)));
            pyramid.triangles.push_back(Triangle(top, Point(pyramid.lowestPoint.x, pyramid.lowestPoint.y + pyramid.width, pyramid.lowestPoint.z),
                                                 Point(pyramid.lowestPoint.x, pyramid.lowestPoint.y, pyramid.lowestPoint.z)));
            pyramid.triangles.push_back(Triangle(top, Point(pyramid.lowestPoint.x + pyramid.width, pyramid.lowestPoint.y, pyramid.lowestPoint.z),
                                                 Point(pyramid.lowestPoint.x + pyramid.width, pyramid.lowestPoint.y + pyramid.width, pyramid.lowestPoint.z)));
            pyramid.triangles.push_back(Triangle(top, Point(pyramid.lowestPoint.x + pyramid.width, pyramid.lowestPoint.y + pyramid.width, pyramid.lowestPoint.z),
                                                  Point(pyramid.lowestPoint.x, pyramid.lowestPoint.y + pyramid.width, pyramid.lowestPoint.z)));
            //pyramid.printPyramid();
            pyramids.push_back(pyramid);
        }
    }

    double lightSrcCnt;
    scene >> lightSrcCnt;
    //cout << lightSrcCnt << endl;
    for(int i = 0; i < lightSrcCnt; i++)
    {
        LightSource ls;
        scene >> ls.position.x >> ls.position.y >> ls.position.z >> ls.falloff;
        //ls.printLightSource();
        lightSources.push_back(ls);
    }

    double spotlightSrcCnt;
    scene >> spotlightSrcCnt;
    //cout << spotlightSrcCnt << endl;
    for(int i = 0; i < spotlightSrcCnt; i++)
    {
        SpotLightSource ls;
        scene >> ls.position.x >> ls.position.y >> ls.position.z >> ls.falloff;
        scene >> ls.lookAt.x >> ls.lookAt.y >> ls.lookAt.z;
        ls.lookAt = ls.lookAt - ls.position;
        normalize(ls.lookAt);
        scene >> ls.cutoffAngle;
        //ls.printLightSource();
        spotLightSources.push_back(ls);
    }
}

void init(){
    drawaxes = 1;
    drawgrid = 0;
	pos.x = -200;
	pos.y = 50;
	pos.z = 50;
	l.x = 1;
	l.y = 0;
	l.z = 0;
	r.x = 0;
	r.y = -1;
	r.z = 0;
	u.x = 0;
	u.y = 0;
	u.z = 1;
	readInput();
	loadTextureImage();

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(fovY, ar, n,	f);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
