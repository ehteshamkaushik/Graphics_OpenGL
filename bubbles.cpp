#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <iostream>

#include <windows.h>
#include <glut.h>
using namespace std;
#define pi (2*acos(0.0))

class Point
{
public:
    double x, y;
    Point()
    {

    }
    Point(double x, double y)
    {
        this->x = x;
        this->y = y;
    }

    void operator = (const Point &p)
    {
        x = p.x;
        y = p.y;
    }

    Point operator+(const Point& p)
    {
        Point point;
        point.x = this->x + p.x;
        point.y = this->y + p.y;
        return point;
    }
    Point operator-(const Point& p)
    {
        Point point;
        point.x = this->x - p.x;
        point.y = this->y - p.y;
        return point;
    }
    Point operator*(const double& m)
    {
        Point point;
        point.x = this->x * m;
        point.y = this->y * m;
        return point;
    }

    Point operator/(const double& m)
    {
        Point point;
        point.x = this->x / m;
        point.y = this->y / m;
        return point;
    }

    Point operator-(const double& m)
    {
        Point point;
        point.x = this->x - m;
        point.y = this->y - m;
        return point;
    }
    Point operator+(const double& m)
    {
        Point point;
        point.x = this->x + m;
        point.y = this->y + m;
        return point;
    }
};

Point pos, l, r, u;

void printPoint(Point p)
{
    cout << p.x << " " << p.y << endl;
}

Point py, vy, pg, vg;

double calcDist(Point p, Point q)
{
    double d = sqrt((p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y));
    return d;
}

bool checkDist(Point p)
{
    double d = sqrt(p.x*p.x + p.y*p.y);
    if (d + 10 > 100)
    {
        return true;
    }
    return false;
}

void normalize(Point &n)
{
    double v = sqrt(n.x*n.x + n.y*n.y);
    n = n / v;
}

void rotateVector(Point& p, double a)
{
    Point perp;
    perp.x = -p.y;
    perp.y = p.x;
    p = p*cos(a*pi/180) + perp*sin(a*pi/180);
    normalize(p);
}



void reflect(Point &p, Point i, Point o)
{
    Point n;
    n = o - i;
    normalize(n);
    double v = p.x*n.x + p.y*n.y;
    v = 2 * v;
    n = n * v;
    p = p - n;
    normalize(p);
}
void drawCircle(double radius,int segments)
{
    int i;
    Point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
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

void drawArrow(Point p1, Point p2)
{
    Point t(p1.x+7*p2.x, p1.y+7*p2.y);
    Point perp;
    perp.x = -p2.y;
    perp.y = p2.x;
    normalize(perp);
    Point t2(t.x+3*perp.x, t.y+3*perp.y);
    Point t3(t.x-3*perp.x, t.y-3*perp.y);
    glBegin(GL_LINES);
    {
        glVertex3f(p1.x,p1.y,0);
        glVertex3f(p1.x + 10*p2.x, p1.y + 10 * p2.y,0);
    }
    glEnd();
    glBegin(GL_TRIANGLES);
    {
        glColor3f(1, 0, 0);
        glVertex3f(p1.x + 10*p2.x, p1.y + 10 * p2.y,0);
        glVertex3f(t3.x, t3.y,0);
        glVertex3f(t2.x, t2.y,0);

    }
    glEnd();
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key

			break;
		case GLUT_KEY_UP:		// up arrow key

			break;

		case GLUT_KEY_RIGHT:
		    rotateVector(vy, -4);
			break;
		case GLUT_KEY_LEFT:
		    rotateVector(vy, 4);
			break;
		case GLUT_KEY_PAGE_UP:
			break;
		case GLUT_KEY_PAGE_DOWN:
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

	gluLookAt(0,0,150,	0,0,0,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

    glColor3f(1,1,1);
    drawCircle(100,98);

    glPushMatrix();
    {
        glColor3f(1,1,1);
        drawArrow(py, vy);
        glTranslatef(py.x,py.y,0);
        glColor3f(1,1,0);
        drawCircle(10,20);
    }
    glPopMatrix();
    glColor3f(1,1,1);
    drawArrow(pg, vg);
    glTranslatef(pg.x,pg.y,0);
    glColor3f(0,1,0);
    drawCircle(10,20);


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//codes for any changes in Models, Camera
	Point val;
    py = py + vy * 0.05;
    if (checkDist(py))
    {
        Point o(0, 0);
        Point i(py.x + 10 * vy.x, py.y+ 10*vy.y);
        reflect(vy, i, o);
    }
    pg = pg + vg * 0.05;
    if (checkDist(pg))
    {
        Point o(0, 0);
        Point i(pg.x + 10 * vg.x, pg.y+ 10*vg.y);
        reflect(vg, i, o);
    }

    if (calcDist(py, pg) <= 20)
    {
        Point i;
        i = py + pg;
        i = i / 2;
        reflect(vy, i, py);
        reflect(vg, i, pg);
    }
	glutPostRedisplay();
}

void init(){
	//codes for initialization
    py.x = 20;
    py.y = 20;
    pg.x = -20;
    pg.y = -20;
    vy.x = 1;
    vy.y = 2;
    normalize(vy);
    vg.x = -2;
    vg.y = 1;
    normalize(vg);


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

