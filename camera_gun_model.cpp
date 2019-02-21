#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>

#include <windows.h>
#include <glut.h>
#include <vector>
using namespace std;
#define pi (2*acos(0.0))

int drawaxes, drawgrid;

double uAngle, rAngle, mAngle, tAngle;

double cnt = 0;

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
};

Point pos, l, r, u;
vector<Point> d;

void printPoint(Point p)
{
    cout << p.x << " " << p.y << " " << p.z << endl;
}

void rotateVector(Point &l, Point &r, Point &u, double a)
{
    Point temp;
    temp = l*cos(a*pi/180) + u*sin(a*pi/180);
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

void drawFire(Point p, double a)
{
    glColor3f(1, 0, 0);
	glBegin(GL_QUADS);{
		glVertex3f(390, -p.y+a, -p.z+a);
		glVertex3f(390, -p.y+a,-p.z-a);
		glVertex3f(390, -p.y-a,-p.z-a);
		glVertex3f(390, -p.y-a, -p.z+a);
	}glEnd();
}


void drawSquare(double a)
{
    glColor3f(0.7,0.7,0.7);
	glBegin(GL_QUADS);{
		glVertex3f(400, a, a);
		glVertex3f(400, a,-a);
		glVertex3f(400, -a,-a);
		glVertex3f(400, -a, a);
	}glEnd();
}



void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}


void drawGunBack(double radius,int segments)
{
    int i, l = 40;
    Point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*pi/2);
        points[i].y=radius*sin(((double)i/(double)segments)*pi/2);
    }

    glTranslatef(-radius, 0, 0);
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }

    glBegin(GL_LINES);
    {
        glVertex3f(0,0+radius,0);
        glVertex3f(0-l,0+radius,0);
    }
    glEnd();

    glPushMatrix();
    {
        glTranslatef(-l, 0, 0);
        glRotatef(90, 0, 0, 1);
        //draw segments using generated points
        for(i=0;i<segments;i++)
        {
            glBegin(GL_LINES);
            {
                glVertex3f(points[i].x,points[i].y,0);
                glVertex3f(points[i+1].x,points[i+1].y,0);
            }
            glEnd();
        }
    }
    glPopMatrix();
}


void drawGunFront(double radius,int segments)
{
    int i, l = 90;
    Point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*pi/2);
        points[i].y=radius*sin(((double)i/(double)segments)*pi/2);
    }

    glTranslatef(-radius, 0, 0);
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }

    glBegin(GL_LINES);
    {
        glVertex3f(0,0+radius,0);
        glVertex3f(0-l,0+radius,0);
    }
    glEnd();

    glPushMatrix();
    {
        glTranslatef(-l, 0, 0);
        glRotatef(90, 0, 0, 1);
        //draw segments using generated points
        for(i=0;i<segments;i++)
        {
            glBegin(GL_LINES);
            {
                glVertex3f(2*radius - points[i].x, points[i].y,0);
                glVertex3f(2*radius - points[i+1].x, points[i+1].y,0);
            }
            glEnd();
        }
    }
    glPopMatrix();
}

void normalize(Point &n)
{
    double v = sqrt(n.x*n.x + n.y*n.y+n.z*n.z);
    n = n / v;
}

Point fireGun(double a)
{
    double xyProj1 = 60 * cos(uAngle*pi/180);
    double zc1 = 60 * sin(uAngle*pi/180);
    double xc1 = xyProj1*cos(rAngle*pi/180);
    double yc1 = xyProj1*sin(rAngle*pi/180);

    double xyProj2 = (160) * cos((uAngle+mAngle)*pi/180);
    double zc2 = (160) * sin((uAngle+mAngle)*pi/180);
    double xc2 = xyProj2*cos(rAngle*pi/180);
    double yc2 = xyProj2*sin(rAngle*pi/180);
    Point p(xc1, yc1, zc1);
    Point v(xc2-xc1, yc2-yc1, zc2-zc1);
    normalize(v);
    double t = (400-p.x)/v.x;
    Point i(p.x+t*v.x, p.y+t*v.y, p.z+t*v.z);
    if(i.y < a && i.y > -a && i.z < a && i.z > -a)
    {
        d.push_back(i);
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
        case 'w':
            if (rAngle <= 45)
                rAngle += 2;
			break;
        case 'q':
            if (rAngle >= -45)
                rAngle -= 2;
			break;
        case 'e':
            if (uAngle >= -45)
                uAngle -= 2;
			break;
        case 'r':
            if (uAngle <= 45)
                uAngle += 2;
			break;

        case 'd':
            if (tAngle <= 45)
                tAngle += 2;
			break;
        case 'f':
            if (tAngle >= -45)
                tAngle -= 2;
			break;
        case 'a':
            if (mAngle >= -45)
                mAngle -= 2;
			break;
        case 's':
            if (mAngle <= 45)
                mAngle += 2;
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
                fireGun(100);
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

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawSquare(100);
	for(int i = 0; i < d.size(); i++)
    {
        drawFire(d[i], 2);
    }
	glRotatef(180, 0, 1, 0);
	glRotatef(rAngle, 0, 0, 1);
	glRotatef(uAngle, 0, 1, 0);

	glPushMatrix();
	{
        glTranslatef(-60, 0, 0);
        glRotatef(mAngle, 0, 1, 0);
        glRotatef(tAngle, 1, 0, 0);
        int color = 1;
        for(int i = 0; i < 360*2; i++)
        {
            if(i % 20 == 0)
            {
                color = 1-color;
                glColor3f(color,color,color);
            }
            glPushMatrix();
            {
                glRotatef(i/2, 1, 0, 0);
                drawGunFront(5,24);
            }
            glPopMatrix();
        }
	}
	glPopMatrix();
	int color = 1;
    for(int i = 0; i < 360*2; i++)
    {
        if(i % 20 == 0)
        {
            color = 1-color;
            glColor3f(color,color,color);
        }
        glPushMatrix();
        {
            glRotatef(i/2, 1, 0, 0);
            drawGunBack(10,24);
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

void init(){
    drawaxes = 1;
    drawgrid = 0;
	pos.x = -100;
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
	uAngle = 0;
	rAngle = 0;
	mAngle = 0;
	tAngle = 0;

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
	gluPerspective(80,	1,	1,	1000.0);
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

