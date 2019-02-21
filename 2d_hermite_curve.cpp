#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<vector>

#include <windows.h>
#include <glut.h>

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

Point cp[200];
int cpidx, mode, u, g, a;
Point up;
double bh[4][4] = {{2, -2, 1, 1},
                    {-3, 3, -2, -1},
                    {0, 0, 1, 0},
                    {1,0, 0, 0}};
double c_x[4] = {0, 0, 0, 0};
double c_y[4] = {0, 0, 0, 0};
double g_x[4] = {0, 0, 0, 0};
double g_y[4] = {0, 0, 0, 0};

int idx, tidx;
double delay;

std::vector<Point> curvePoints;

void getCoeff()
{
    for(int i = 0; i < 4; i++)
    {
        double sum_x = 0, sum_y = 0;
        for(int j = 0; j < 4; j++)
        {
            sum_x += (bh[i][j] * g_x[j]);
            sum_y += (bh[i][j] * g_y[j]);
        }
        c_x[i] = sum_x;
        c_y[i] = sum_y;
    }
}

void drawSquare()
{
    glBegin(GL_QUADS);
    {
        glVertex3d( 3,  3, 0);
        glVertex3d( 3, -3, 0);
        glVertex3d(-3, -3, 0);
        glVertex3d(-3,  3, 0);
    }
    glEnd();
}

void drawCircle(double radius,int segments)
{
    int i;
    Point points[100];
    glColor3f(0, 1, 1);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,0);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

double calcDist(Point p, Point q)
{
    double d = sqrt((p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y));
    return d;
}


void normalize(Point &n)
{
    double v = sqrt(n.x*n.x + n.y*n.y);
    n = n / v;
}

void drawArrow(Point p1, Point p2, Point p3)
{
    double d;
    d = calcDist(p1, p3);
    d = d*20/100;
    Point t(p1.x-d*p2.x, p1.y-d*p2.y);
    Point perp;
    perp.x = -p2.y;
    perp.y = p2.x;
    normalize(perp);
    Point t2(t.x+7*perp.x, t.y+7*perp.y);
    Point t3(t.x-7*perp.x, t.y-7*perp.y);

    glBegin(GL_TRIANGLES);
    {
        glColor3f(1, 0, 0);
        glVertex3f(p1.x-3*p2.x , p1.y-3*p2.y,0);
        glVertex3f(t3.x, t3.y,0);
        glVertex3f(t2.x, t2.y,0);
    }
    glEnd();
}

void drawLine(Point p1, Point p2)
{
    Point v;
    v.x = p2.x - p1.x;
    v.y = p2.y - p1.y;
    normalize(v);

    drawArrow(p2, v, p1);
    glColor3f(1, 1, 1);
    glBegin(GL_LINES);
    {
        glVertex3d( p1.x,  p1.y, 0);
        glVertex3d( p2.x, p2.y, 0);
    }
    glEnd();

}

void generateCurvePoints(Point p1, Point p2, Point p4, Point p3)
{
    Point r1(p2.x - p1.x, p2.y - p1.y);
    Point r4(p3.x - p4.x, p3.y - p4.y);
    g_x[0] = p1.x;
    g_x[1] = p4.x;
    g_x[2] = r1.x;
    g_x[3] = r4.x;
    g_y[0] = p1.y;
    g_y[1] = p4.y;
    g_y[2] = r1.y;
    g_y[3] = r4.y;
    getCoeff();
    int n = 100;
    double d = 1.0/n*1.0;
    double delay = 0;
    double f_x = c_x[3];
    double f_y = c_y[3];
    double del_f_x = c_x[0]*d*d*d + c_x[1]*d*d + c_x[2]*d;
    double del_f_y = c_y[0]*d*d*d + c_y[1]*d*d + c_y[2]*d;
    double del_f2_x = 6 * c_x[0]*d*d*d + 2 * c_x[1]*d*d;
    double del_f2_y = 6 * c_y[0]*d*d*d + 2 * c_y[1]*d*d;
    double del_f3_x = 6 * c_x[0]*d*d*d;
    double del_f3_y = 6 * c_y[0]*d*d*d;
    curvePoints.push_back(Point(f_x, f_y));
    for(int i = 1; i <= n; i++)
    {
        f_x += del_f_x;
        del_f_x += del_f2_x;
        del_f2_x += del_f3_x;
        f_y += del_f_y;
        del_f_y += del_f2_y;
        del_f2_y += del_f3_y;
        curvePoints.push_back(Point(f_x, f_y));
    }
}

void drawCurve()
{
    int n = curvePoints.size();
    glColor3f(1, 1, 1);
    for(int i = 0; i < n - 1; i++)
    {
        glBegin(GL_LINES);
        {
            glVertex3d( curvePoints[i].x,  curvePoints[i].y, 0);
            glVertex3d( curvePoints[i+1].x, curvePoints[i+1].y, 0);
        }
        glEnd();
    }
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

        case 'g':
            g = 1 - g;
			break;
        case 'u':
            if(u == 0 && mode >= 2)
            {
                mode = 3;
                delay = 1;
                tidx = 0;
                a = 0;
            }
            else if(u == 1 && mode == 3)
            {
                mode = 2;
                u = 1 - u;
            }
			break;
        case 'a':
            if(a == 0 && mode == 2)
            {
                mode = 4;
                a = 1 - a;

            }
            else if(a == 1 && mode == 4)
            {
                delay = 1;
                tidx = 0;
                mode = 2;
                a = 1 - a;
            }
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
			break;
		case GLUT_KEY_LEFT:
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

int update_state;
void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){
                if(mode == 0)
                {
                    cp[cpidx].x = (double)x;
                    cp[cpidx].y = (double)(600 - y);
                    cpidx++;
                }
                if(mode == 3)
                {
                    if(update_state == 0)
                    {
                        up.x = (double)x;
                        up.y = (double)(600 - y);
                        update_state = 1;
                    }
                    else if(update_state == 1)
                    {
                        up.x = (double)x;
                        up.y = (double)(600 - y);
                        update_state = 2;
                    }
                }
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			if(state == GLUT_DOWN){
                if(cpidx > 3 && cpidx % 2 == 0)
                mode = 1;
			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void drawMain()
{
    int i;
    int r = 0;

    for (i = 0; i < cpidx; i++)
    {
        glColor3f(r, 1, 0);
        glPushMatrix();
        {
            r = 1 - r;
            glTranslatef(cp[i].x, cp[i].y, 0);
            drawSquare();
        }
        glPopMatrix();
        if(i % 2 == 1)
        {
            glColor3f(1, 1, 1);
            drawLine(cp[i-1], cp[i]);
        }
    }
}

void curveInit()
{
    curvePoints.clear();
    if(cpidx > 3 && cpidx % 2 == 0)
    {
        for (int j = 0; j < cpidx; j += 2)
        {
            generateCurvePoints(cp[j%(cpidx)], cp[(j+1)%(cpidx)], cp[(j+2)%(cpidx)], cp[(j+3)%(cpidx)]);
        }
    }
}

void pathTraveller()
{
    int n = curvePoints.size();
    glPushMatrix();
    {
        glTranslatef(curvePoints[tidx % n].x, curvePoints[tidx % n].y, 0);
        drawCircle(5, 10);
    }
    glPopMatrix();
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
	//gluLookAt(150*cos(cameraAngle), 150*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(0,0,0,	0,0,-1,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	if(mode == 0)
    {
        drawMain();
    }


    if (mode == 1)
    {
        curveInit();
        mode = 2;
    }
    if (mode == 2)
    {
        if(g == 0)
        {
            drawMain();
        }
        drawCurve();
    }

    if (mode == 3)
    {
        drawMain();
        drawCurve();
        if (update_state == 1)
        {
            idx = -1;
            double min_dist = 10000000;
            for(int i = 0; i < cpidx; i++)
            {
                double curr_dist = calcDist(up, cp[i]);
                if(curr_dist < min_dist)
                {
                    idx = i;
                    min_dist = curr_dist;
                }
            }

            glPushMatrix();
            {
                glTranslatef(cp[idx].x, cp[idx].y, 0);
                drawCircle(10, 10);
            }
            glPopMatrix();
        }
        if(update_state == 2)
        {
            cp[idx].x = up.x;
            cp[idx].y = up.y;
            curveInit();

            u = 0;
            g = 0;
            update_state = 0;
            mode = 1;
        }
    }

    if(mode == 4)
    {
        if(g == 0)
        {
            drawMain();
        }
        drawCurve();
        pathTraveller();
    }


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
    if (mode == 4)
    {
        delay += 1;

        if((int)delay % 10 == 0)
        {
            delay = 1;
            tidx += 1;
        }

    }

	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization

	cpidx = 0;
    mode = 0;
    u = 0;
    g = 0;
    a = 0;
    update_state = 0;
    idx = 0;
    delay = 1;
    tidx = 0;
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
	gluOrtho2D(0, 800, 0, 600);
	//gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(800, 600);
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
