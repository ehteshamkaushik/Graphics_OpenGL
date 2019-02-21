#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include <algorithm>
#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))
#define epsilon (1.0e-6)

class homogeneous_point
{
public:
    double x, y, z, w;

    // set the three coordinates, set w to 1
    homogeneous_point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }

    /*
    default constructor. does nothing. allows declarations like below:
        matrix m;
    therefore, usage is dangerous
    */
    homogeneous_point() {
    }

    // constructs a homogeneous point with given coordinates. forces w to be 1.0
    // if w is zero, raises error
    homogeneous_point(double x, double y, double z, double w)
    {
        assert (w != 0);
        this->x = x/w;
        this->y = y/w;
        this->z = z/w;
        this->w = 1;
    }

    // adds two points. returns a point forcing w to be 1.0
    homogeneous_point operator+ (const homogeneous_point& point)
    {
        double x = this->x + point.x;
        double y = this->y + point.y;
        double z = this->z + point.z;
        double w = this->w + point.w;
        homogeneous_point p(x, y, z, w);
        return p;
    }

    // subtracts one point from another. returns a point forcing w to be 1.0
    homogeneous_point operator- (const homogeneous_point& point)
    {
        double x = this->x - point.x;
        double y = this->y - point.y;
        double z = this->z - point.z;
        double w = this->w - point.w;
        homogeneous_point p(x, y, z, w);
    }

    // Print the coordinates of a point. exists for testing purpose.
    void print()
    {
        cout << "Point: " << endl;
        cout << x << " " << y << " " << z << " " << w << endl;
    }

};


class Vector
{
public:
    double x, y, z;

    // constructs a vector with given components
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // keeps the direction same. recalculates the vector to be unit.
    void normalize()
    {
        double r = sqrt(x*x + y*y + z*z);
        x = x / r;
        y = y / r;
        z = z / r;
    }

    // add two vectors
    Vector operator+(const Vector& v)
    {
        Vector v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }

    // subtract one vector from another
    Vector operator-(const Vector& v)
    {
        Vector v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

    // scale a vector with a given coefficient
    Vector operator* (double m)
    {
        Vector v(x*m, y*m, z*m);
        return v;
    }

    // get the dot product of two vectors
    static double dot(Vector a, Vector b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    // get the cross product of two vectors
    static Vector cross(Vector a, Vector b)
    {
        Vector v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
        return v;
    }

    // print a vector. only for testing purposes.
    void print ()
    {
        cout << "Vector" << endl;
        cout << x << " " << y << " " << z << endl;
    }
};


/*
The matrices are forced to be 4x4. This is because in this assignment, we will deal with points in triangles.
Maximum # of points that we will deal with at once is 3. And all the standard matrices are 4x4 (i.e. scale, translation, rotation etc.)
*/
class matrix
{
public:
    double values[4][4];
    int num_rows, num_cols;

    // only set the number of rows and cols
    matrix(int rows, int cols)
    {
        assert (rows <= 4 && cols <= 4);
        num_rows = rows;
        num_cols = cols;
    }

    // prepare an nxn square matrix
    matrix(int n)
    {
        assert (n <= 4);
        num_rows = num_cols = n;
    }

    // prepare and return an identity matrix of size nxn
    static matrix make_identity(int n)
    {
        assert (n <= 4);
        matrix m(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    m.values[i][j] = 1;
                else
                    m.values[i][j] = 0;
            }
        }
        return m;
    }

    // print the matrix. exists for testing purposes
    void print()
    {
        cout << "Matrix:" << endl;
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                cout << values[i][j] << "\t";
            }
            cout << endl;
        }
    }

    // add the two matrices. Raise error if dimension mismatches
    matrix operator+ (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] + m.values[i][j];
            }
        }
        return m1;
    }

    // subtract a matrix from another. raise error if dimension mismatches
    matrix operator- (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] - m.values[i][j];
            }
        }
        return m1;
    }

    // multiply two matrices. allows statements like m1 = m2 * m3; raises error is dimension mismatches
    matrix operator* (const matrix& m)
    {
        assert (this->num_cols == m.num_rows);
        matrix m1(this->num_rows, m.num_cols);

        for (int i = 0; i < m1.num_rows; i++) {
            for (int j = 0; j < m1.num_cols; j++) {
                double val = 0;
                for (int k = 0; k < this->num_cols; k++) {
                    val += this->values[i][k] * m.values[k][j];
                }
                m1.values[i][j] = val;
            }
        }
        return m1;
    }

    // multiply a matrix with a constant
    matrix operator* (double m)
    {
        matrix m1(this->num_rows, this->num_cols);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m1.values[i][j] = m * this->values[i][j];
            }
        }
        return m1;
    }

    // multiply a 4x4 matrix with a homogeneous point and return the resulting point.
    // usage: homogeneous_point p = m * p1;
    // here, m is a 4x4 matrix, intended to be the transformation matrix
    // p1 is the point on which the transformation is being made
    // p is the resulting homogeneous point
    homogeneous_point operator* (const homogeneous_point& p)
    {
        assert (this->num_rows == this->num_cols && this->num_rows == 4);

        matrix m(4, 1);
        m.values[0][0] = p.x;
        m.values[1][0] = p.y;
        m.values[2][0] = p.z;
        m.values[3][0] = p.w;

        matrix m1 = (*this)*m;
        homogeneous_point p1(m1.values[0][0], m1.values[1][0], m1.values[2][0], m1.values[3][0]);
        return p1;
    }

    // return the transpose of a matrix
    matrix transpose()
    {
        matrix m(num_cols, num_rows);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m.values[j][i] = values[i][j];
            }
        }
        return m;
    }

};

/*
A simple class to hold the color components, r, g, b of a certain shade.
*/
class color {
public:
    double r, g, b;
    color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    color() {
    }
};


double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fov_x, fov_y, aspectRatio, near, far;
color backgroud;
int screen_x, screen_y;
int triangle_count = 0;
vector <color> initialColor;
vector <color> clippedColor;



void scan_convert() {
    ifstream stage3;
    stage3.open("stage3.txt");
    cout << "Stage 4\n";

    color** pixels = new color*[screen_x];
    double** zs = new double*[screen_x];
    for (int i = 0; i < screen_x; i++) {
        pixels[i] = new color [screen_y];
        for (int j = 0; j < screen_y; j++) {
            pixels[i][j] = backgroud;
        }
        zs[i] = new double [screen_y];
        for (int j = 0; j < screen_y; j++) {
            zs[i][j] = +20; // a very large value intended as +INFINITY
        }
    }

    // perform scan conversion, populate the 2D array pixels
    // the array zs is the z-buffer.

    for(int i = 0; i < triangle_count; i++)
    {
        //cout << i << endl;
        vector<homogeneous_point> triangle;
        for(int j = 0; j < 3; j++)
        {
            double x, y, z;
            stage3 >> x >> y >> z;
            homogeneous_point p(x, y, z);
            triangle.push_back(p);
        }
        sort(triangle.begin(), triangle.end(), [](const homogeneous_point& lhs, const homogeneous_point& rhs)
             {
                 return lhs.y > rhs.y;
             });
        homogeneous_point p1 = triangle.at(0);
        homogeneous_point p2 = triangle.at(1);
        homogeneous_point p3 = triangle.at(2);
        double diff_y = 2.0/screen_y*1.0;
        double diff_x = 2.0/screen_x*1.0;
        double ys = p1.y;

        while(true)
        {
            if(ys > 1)
                continue;
            else if(ys < -1)
                break;
            double za, zb, xa, xb, zp;
            if(ys > p2.y)
            {
                za = p1.z - (p1.z - p2.z)*(p1.y - ys)/(p1.y - p2.y);
                xa = p1.x - (p1.x - p2.x)*(p1.y - ys)/(p1.y - p2.y);
            }
            else
            {
                za = p2.z - (p2.z - p3.z)*(p2.y - ys)/(p2.y - p3.y);
                xa = p2.x - (p2.x - p3.x)*(p2.y - ys)/(p2.y - p3.y);
            }
            zb = p1.z - (p1.z - p3.z)*(p1.y - ys)/(p1.y - p3.y);
            xb = p1.x - (p1.x - p3.x)*(p1.y - ys)/(p1.y - p3.y);
            if(xa > xb)
            {
                swap(xa, xb);
                swap(za, zb);
            }
            //cout << "ZA : " << za << " ZB : " << zb << endl;
            //cout << "XA : " << xa << " XB : " << xb << endl;
            for (double xp = xa; xp < xb; xp+=diff_x)
            {
                if(xp < -1)
                    continue;
                else if(xp > 1)
                    break;
                zp = zb - (zb - za)*(xb - xp)/(xb-xa);
                int pixel_y = (int)((1+ys)/diff_y);
                int pixel_x = (int)((1+xp)/diff_x);
                if(zp+epsilon < zs[pixel_x][screen_y - pixel_y])
                {
                    zs[pixel_x][screen_y - pixel_y] = zp;
                    pixels[pixel_x][screen_y - pixel_y] = clippedColor.at(i);
                }
            }
            ys = ys - diff_y;
            if(ys <= p3.y)
                break;
        }
    }


    // the following code generates a bmp image. do not change this.
    bitmap_image image(screen_x, screen_y);
    for (int x = 0; x < screen_x; x++) {
        for (int y = 0; y < screen_y; y++) {
            image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
        }
    }
    image.save_image("out.bmp");

    // free the dynamically allocated memory

}


void stage3()
{
    if (near == far) return;
    ifstream stage2;
    ofstream stage3;
    stage2.open ("stage2.txt");
    stage3.open ("stage3.txt");
    stage3 << std::fixed;
    stage3 << std::setprecision(7);

    // process input from stage2 and write to stage3
    fov_x = fov_y * aspectRatio;
    double t = near * tan((fov_y/2) * (pi/180));
    double r = near * tan((fov_x/2) * (pi/180));

    matrix P = matrix(4);

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            P.values[i][j] = 0;
        }
    }
    P.values[0][0] = near/r;

    P.values[1][1] = near/t;

    P.values[2][2] = -(far + near)/(far - near);
    P.values[2][3] = -(2 * far * near)/(far - near);

    P.values[3][2] = -1;

    vector<homogeneous_point> clippedPoints;

    for(int i = 0; i < triangle_count; i++)
    {
        cout << "Triangle : " << i << endl;
        vector<homogeneous_point> inputlist;
        vector<homogeneous_point> outputlist;
        for(int j = 0; j < 3; j++)
        {
            double x, y, z;
            stage2 >> x >> y >> z;
            homogeneous_point p(x, y, z);
            inputlist.push_back(p);
        }
        double edge = -near;
        double j = 0;
        double n = inputlist.size() - 1;
        homogeneous_point I = inputlist.at(n);
        homogeneous_point F = inputlist.at(0);

        while(true)
        {
            if (F.z <= edge && I.z <= edge)
            {
                outputlist.push_back(F);
            }

            else if(I.z <= edge && F.z > edge)
            {
                homogeneous_point intersect;
                intersect.z = edge;
                intersect.x = I.x - (I.x - F.x)*(I.z - edge)/(I.z - F.z);
                intersect.y = I.y - (I.y - F.y)*(I.z - edge)/(I.z - F.z);
                outputlist.push_back(intersect);
            }
            else if(I.z > edge && F.z <= edge)
            {
                homogeneous_point intersect;
                intersect.z = edge;
                intersect.x = I.x - (I.x - F.x)*(I.z - edge)/(I.z - F.z);
                intersect.y = I.y - (I.y - F.y)*(I.z - edge)/(I.z - F.z);
                outputlist.push_back(intersect);
                outputlist.push_back(F);
            }
            j++;
            if(j > n)
                break;
            I = inputlist.at(j-1);
            F = inputlist.at(j);
        }

        if(outputlist.size() == 0)
            continue;
        cout << outputlist.size() << endl;
        inputlist.clear();
        for (int j = 0; j < outputlist.size(); j++)
        {
            inputlist.push_back(outputlist.at(j));
        }
        outputlist.clear();

        edge = -far;
        j = 0;
        n = inputlist.size() - 1;
        I = inputlist.at(n);
        F = inputlist.at(0);

        while(true)
        {
            if (F.z >= edge && I.z >= edge)
            {
                outputlist.push_back(F);
            }

            else if(I.z >= edge && F.z < edge)
            {
                homogeneous_point intersect;
                intersect.z = edge;
                intersect.x = I.x - (I.x - F.x)*(I.z - edge)/(I.z - F.z);
                intersect.y = I.y - (I.y - F.y)*(I.z - edge)/(I.z - F.z);
                outputlist.push_back(intersect);
            }
            else if(I.z < edge && F.z >= edge)
            {
                homogeneous_point intersect;
                intersect.z = edge;
                intersect.x = I.x - (I.x - F.x)*(I.z - edge)/(I.z - F.z);
                intersect.y = I.y - (I.y - F.y)*(I.z - edge)/(I.z - F.z);
                outputlist.push_back(intersect);
                outputlist.push_back(F);
            }
            j++;
            if(j > n)
                break;
            I = inputlist.at(j-1);
            F = inputlist.at(j);
        }
        if(outputlist.size() == 0)
            continue;
        for (int j = 1; j <= outputlist.size() - 2; j++)
        {
            clippedPoints.push_back(outputlist.at(0));
            clippedPoints.push_back(outputlist.at(j));
            clippedPoints.push_back(outputlist.at(j+1));
            clippedColor.push_back(initialColor.at(i));
        }
    }

    for (int i = 0; i < clippedPoints.size(); i++)
    {
        homogeneous_point p = clippedPoints.at(i);
        homogeneous_point pNew = P * p;
        stage3 << pNew.x << " " << pNew.y << " " << pNew.z << '\n';
        if((i+1) % 3 == 0)
            stage3 << '\n';
    }
    triangle_count = clippedPoints.size()/3;
    stage3.close();
    stage2.close();
}

void stage2()
{
    ifstream stage1;
    ofstream stage2;
    stage1.open ("stage1.txt");
    stage2.open ("stage2.txt");
    stage2 << std::fixed;
    stage2 << std::setprecision(7);

    // collect input from stage1 and process, write output to stage2

    Vector eye(eye_x, eye_y, eye_z);
    Vector look(look_x, look_y, look_z);
    Vector up(up_x, up_y, up_z);

    Vector l = look - eye;
    l.normalize();
    Vector r = Vector::cross(l, up);
    r.normalize();
    Vector u = Vector::cross(r, l);
    u.normalize();

    matrix T = matrix::make_identity(4);
    T.values[0][3] = -eye_x;
    T.values[1][3] = -eye_y;
    T.values[2][3] = -eye_z;

    matrix R = matrix::make_identity(4);
    R.values[0][0] = r.x;
    R.values[1][0] = u.x;
    R.values[2][0] = -l.x;

    R.values[0][1] = r.y;
    R.values[1][1] = u.y;
    R.values[2][1] = -l.y;

    R.values[0][2] = r.z;
    R.values[1][2] = u.z;
    R.values[2][2] = -l.z;

    matrix V = R * T;

    for(int i = 0; i < triangle_count; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            double x, y, z;
            stage1 >> x >> y >> z;
            homogeneous_point p(x, y, z);
            homogeneous_point pNew = V * p;
            stage2 << pNew.x << " " << pNew.y << " " << pNew.z << '\n';
        }
        stage2 << '\n';
    }

    stage1.close();
    stage2.close();

}

void showstack(stack <matrix> s)
{
    while (!s.empty())
    {
        s.top().print();
        s.pop();
    }
    cout << '\n';
}

Vector R(Vector x, Vector a, double theta)
{
    double redianAngle = theta*pi/180;
    return (x * cos(redianAngle)) + (a * (Vector::dot(a, x) * (1 - cos(redianAngle)))) + (Vector::cross(a, x) * sin(redianAngle));
}

void stage1()
{
    ifstream scene;
    ofstream stage1;
    scene.open ("scene.txt");
    stage1.open ("stage1.txt");
    stage1 << std::fixed;
    stage1 << std::setprecision(7);

    string command = "";

    scene >> eye_x >> eye_y >> eye_z;
    scene >> look_x >> look_y >> look_z;
    scene >> up_x >> up_y >> up_z;
    scene >> fov_y >> aspectRatio >> near >> far;
    scene >> screen_x >> screen_y;
    scene >> backgroud.r >> backgroud.g >> backgroud.b;

    matrix I = matrix::make_identity(4);
    stack<matrix> S;
    stack<int> track;
    S.push(I);


    while(true)
    {
        scene >> command;
        if(command == "triangle")
        {
            triangle_count++;
            matrix transformMatrix = S.top();
            for(int i = 0; i < 3; i++)
            {
                double x, y, z;
                scene >> x >> y >> z;
                homogeneous_point p(x, y, z);
                homogeneous_point pNew = transformMatrix * p;
                stage1 << pNew.x << " " << pNew.y << " " << pNew.z << '\n';
            }
            stage1 << '\n';
            double r, g, b;
            scene >> r >> g >> b;
            color triangleColor(r, g, b);
            initialColor.push_back(triangleColor);
        }
        else if(command == "translate")
        {
            matrix transformMatrix = S.top();
            double tx, ty, tz;
            scene >> tx >> ty >> tz;
            matrix translateMatrix = matrix::make_identity(4);
            translateMatrix.values[0][3] = tx;
            translateMatrix.values[1][3] = ty;
            translateMatrix.values[2][3] = tz;
            matrix newMatrix = transformMatrix * translateMatrix;
            S.push(newMatrix);
        }
        else if(command == "scale")
        {
            matrix transformMatrix = S.top();
            double sx, sy, sz;
            scene >> sx >> sy >> sz;
            matrix scaleMatrix = matrix::make_identity(4);
            scaleMatrix.values[0][0] = sx;
            scaleMatrix.values[1][1] = sy;
            scaleMatrix.values[2][2] = sz;
            matrix newMatrix = transformMatrix * scaleMatrix;
            S.push(newMatrix);
        }
        else if(command == "rotate")
        {
            matrix transformMatrix = S.top();
            double angle, ax, ay, az;
            scene >> angle >> ax >> ay >> az;
            Vector a(ax, ay, az);
            a.normalize();
            Vector i(1, 0, 0);
            Vector j(0, 1, 0);
            Vector k(0, 0, 1);

            Vector c1 = R(i, a, angle);
            Vector c2 = R(j, a, angle);
            Vector c3 = R(k, a, angle);

            matrix rotationMatrix = matrix::make_identity(4);
            rotationMatrix.values[0][0] = c1.x;
            rotationMatrix.values[1][0] = c1.y;
            rotationMatrix.values[2][0] = c1.z;

            rotationMatrix.values[0][1] = c2.x;
            rotationMatrix.values[1][1] = c2.y;
            rotationMatrix.values[2][1] = c2.z;

            rotationMatrix.values[0][2] = c3.x;
            rotationMatrix.values[1][2] = c3.y;
            rotationMatrix.values[2][2] = c3.z;

            matrix newMatrix = transformMatrix * rotationMatrix;
            S.push(newMatrix);

        }
        else if(command == "push")
        {
            track.push(S.size());
        }
        else if(command == "pop")
        {
            int popCount = S.size() - track.top();
            track.pop();
            for (int i = 0; i < popCount; i++)
            {
                S.pop();
            }
        }
        else if(command == "end")
        {
            break;
        }
    }


    // take other commands as input from scene in a loop
    // process accordingly
    // write to stage1

    scene.close();
    stage1.close();

}

int main()
{
    cout << std::fixed;
    cout << std::setprecision(4);

    stage1();
    stage2();
    stage3();
    scan_convert();

    return 0;
}

