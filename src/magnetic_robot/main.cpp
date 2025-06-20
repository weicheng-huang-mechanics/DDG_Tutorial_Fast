/**
 * simDiscretePlate
 * simDiscretePlate stands for "[sim]plified Discrete Plate"
 * NOV 2019
 * Author: Weicheng Huang 
 * */

//This line is for mac
//#include <GLUT/glut.h>

//This is for linux
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <string>
#include "eigenIncludes.h"

// Rod and stepper are included in the world
#include "world.h"
#include "setInput.h"

world myWorld;
int NPTS, NRODS;
ofstream outfile;

static void Key(unsigned char key, int x, int y)
{
  switch (key) // ESCAPE to quit
  {
	case 27:
		exit(0);
  }
}

/* Initialize OpenGL Graphics */
void initGL() 
{
	glClearColor(0.7f, 0.7f, 0.7f, 0.0f); // Set background color to black and opaque
	glClearDepth(10.0f);                   // Set background depth to farthest
	//glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	//glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable smooth shading
	//glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections

	glLoadIdentity();
	gluLookAt(0.05, 0.05, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
	glPushMatrix();

	//glMatrixMode(GL_MODELVIEW);
}

void display(void)
{
	while ( myWorld.simulationRunning() > 0)
	{
		//  Clear screen and Z-buffer
		glClear(GL_COLOR_BUFFER_BIT);

		// draw axis
		double axisLen = 1;
		glLineWidth(0.5);
		
		glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(axisLen, 0.0, 0.0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, axisLen, 0.0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, axisLen);
		glEnd();

		glColor3f(0.1, 0.1, 0.1);
		glLineWidth(3.0);
		
		glBegin(GL_LINES);
		NPTS = myWorld.numStretchingPair();
		for (int i = 0; i < NPTS; i++)
		{
			Vector3d xCurrent_1 = myWorld.getScaledCoordinatePlate(i, 0);
			Vector3d xCurrent_2 = myWorld.getScaledCoordinatePlate(i, 1);
			glVertex3f(xCurrent_1(0), xCurrent_1(1)-0.5, xCurrent_1(2));
			glVertex3f(xCurrent_2(0), xCurrent_2(1)-0.5, xCurrent_2(2));
		}
		glEnd();

		//draw a line
		glColor3f(0.0, 1.0, 0.0);
		glLineWidth(2.0);
		
		for (int n=0; n < NRODS; n++)
		{
			glBegin(GL_LINES);
			NPTS = myWorld.numPoints(n);
			for (int i=0; i < NPTS-1; i++)
			{
				glVertex3f( myWorld.getScaledCoordinate(n, 4*i), myWorld.getScaledCoordinate(n, 4*i+1)-0.5, myWorld.getScaledCoordinate(n, 4*i+2) );
				glVertex3f( myWorld.getScaledCoordinate(n, 4*(i+1)), myWorld.getScaledCoordinate(n, 4*(i+1)+1)-0.5, myWorld.getScaledCoordinate(n, 4*(i+1)+2));
			}
			glEnd();
		}
		
		/*

		glColor3f(0.0, 0.0, 1.0);
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for (int n=0; n < 3 * NRODS; n++)
		{
			Vector3d xCurrent_2 = myWorld.getScaledCoordinateCouple(n, 1);
			glVertex3f(xCurrent_2(0), xCurrent_2(1), xCurrent_2(2));

		}

		glEnd();

		glColor3f(0.0, 0.0, 1.0);
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for (int n=0; n < 3 * NRODS; n++)
		{
			Vector3d xCurrent_1 = myWorld.getScaledCoordinateCouple(n, 0);
			glVertex3f(xCurrent_1(0), xCurrent_1(1), xCurrent_1(2));
		}

		glEnd();

		glColor3f(1.0, 0.0, 0.0);
		glPointSize(2.0);
		glBegin(GL_POINTS);

		Vector3d xCurrent = myWorld.getMidPoint();
		glVertex3f(xCurrent(0), xCurrent(1), xCurrent(2));

		glEnd();

		*/
		
		glFlush();
		
		// Update step
		myWorld.updateTimeStep();
		myWorld.CoutData(outfile);
	}
	exit(1);
}

int main(int argc,char *argv[])
{
	setInput inputData;
	inputData = setInput();
	inputData.LoadOptions(argv[1]);
	inputData.LoadOptions(argc,argv);
	//read input parameters from txt file and cmd

	myWorld = world(inputData);
	myWorld.setPlateStepper();

	myWorld.OpenFile(outfile);

	bool render = myWorld.isRender();

	NPTS = myWorld.numStretchingPair();

	if (render) // if OpenGL visualization is on
	{
		NRODS = myWorld.getnumRods();

		glutInit(&argc,argv);
		glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize (1000, 1000);
		glutInitWindowPosition (100, 100);
		glutCreateWindow ("simDiscretePlate");
		initGL();
		glutKeyboardFunc(Key);
		glutDisplayFunc(display);
		glutMainLoop();
	}	
	else
	{
		while ( myWorld.simulationRunning() > 0)
		{
			myWorld.updateTimeStep(); // update time step
			myWorld.CoutData(outfile); // write data to file
		}
	}

	// Close (if necessary) the data file
	myWorld.CloseFile(outfile);
	
	return 0;
}

