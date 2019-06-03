// vim: tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab

/** @file housePerspective.c
 *
 *  Plots four views of a house using different projection types.
 *  Non-commercial use only.
 *
 *  @author Paulo Roma Cavalcanti
 *  @date 10/01/2017.
 *  @see http://www.glprogramming.com/red/chapter03.html
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<GL/glut.h>

int xsize = 256, ysize = 256;
double ratio = 1.0;

void renderViewportLines (void) {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
 	gluOrtho2D (-1, 1, -1, 1);
	glViewport ( 0, 0, xsize, ysize );

	glBegin (GL_LINES);
		glVertex2f (-1.0, 0.0);
		glVertex2f (1.0, 0.0);

		glVertex2f (0.0, 1.0);
		glVertex2f (0.0, -1.0);
	glEnd(); 
}

void renderHouseFront(double z) {
	glBegin (GL_LINE_LOOP);
		glVertex3f (0.0, 0.0, z);
		glVertex3f (16.0, 0.0, z);
		glVertex3f (16.0, 10.0, z);
		glVertex3f (8.0, 16.0, z);
		glVertex3f (0.0, 10.0, z);
		glVertex3f (0.0, 0.0, z);
	glEnd();
}

void renderScene(void) {
	renderHouseFront(30.0);
	renderHouseFront(54.0);

	glBegin (GL_LINES);
		glColor3f ( 1.0, 0.0, 0.0 );
		glVertex3f (0.0, 0.0, 30.0);
		glVertex3f (0.0, 0.0, 54.0);

		glVertex3f (16.0, 0.0, 30.0);
		glVertex3f (16.0, 0.0, 54.0);

		glColor3f ( 0.0, 0.0, 1.0 );
		glVertex3f (16.0, 10.0, 30.0);
		glVertex3f (16.0, 10.0, 54.0);

		glVertex3f (0.0, 10.0, 30.0);
		glVertex3f (0.0, 10.0, 54.0);

		glColor3f ( 1.0, 1.0, 1.0 );
		glVertex3f (8.0, 16.0, 30.0);
		glVertex3f (8.0, 16.0, 54.0);
	glEnd();
}

/**
 *  Returns the field of view, given window height and
 *  distance to the view point.
 *
 *  @param size window height.
 *  @param distance view point distance.
 *  @return field of view.
 */
double calculateAngle(double size, double distance) {
	double radtheta, degtheta; 

	radtheta = 2.0 * atan2 (size/2.0, distance);
	degtheta = (180.0 * radtheta) / M_PI; // defined in math.h
	return (degtheta);
}

void setLookAt(GLdouble eyex, GLdouble eyey, GLdouble eyez, 
               GLdouble aimx, GLdouble aimy, GLdouble aimz) {
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyex, eyey, eyez, aimx, aimy, aimz, 0.0, 1.0, 0.0);
}

void setFrustum (GLdouble left, GLdouble right, GLdouble bottom, 
                 GLdouble top, GLdouble _near, GLdouble _far) {
	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// The camera is at (0,0) for the frustum.
	// I.e. The center of projection is treated as the center for the frustum.
	// Set the correct perspective.

#ifndef __KEEP_ASPECT__
	glFrustum (left, right, bottom, top, _near, _far);
#else
	double fov = calculateAngle(top-bottom, _near);
	gluPerspective(fov, ratio, _near, _far);
	// printf ("fov = %f\n", fov);
#endif
}

void redraw () {

	int dxsize = 0.5*xsize, dysize = 0.5*ysize;

	// Clear screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	renderViewportLines();

	// Front view: quadrant 1.
 	glViewport (0, dysize, dxsize, dysize );
	setLookAt (8.0,6.0,84.0, 
#ifndef __KEEP_ASPECT__
		8.0,6.0,30.0);
#else
		8.0,8.0,42.0);  // center of the house
#endif
	setFrustum (-9.0, 9.0, -7.0, 11.0, 29, 150);
 	renderScene();

	// Side view: quadrant 2.
 	glViewport ( dxsize, dysize, dxsize, dysize );
	setLookAt (26.0,8.0,42.0, 
#ifndef __KEEP_ASPECT__
		16.0,8.0,42.0);
#else
		8.0,8.0,42.0);  // center of the house
#endif
	setFrustum (-8.0, 8.0, -8.0, 8.0, 5, 35);
 	renderScene();

	// Parallel projection - side view: quadrant 3.
 	glViewport ( 0, 0, dxsize, dysize );
	setLookAt (26.0,8.0,42.0, 
#ifndef __KEEP_ASPECT__
		16.0,8.0,42.0);
	double w = 15.0; 
#else
		8.0,8.0,42.0);  // center of the house
	double w = ratio*15.0; 
#endif
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho (-w, w, -15.0, 15.0, -10, 150);
 	renderScene();

	// 3/4 view: quadrant 4.
 	glViewport ( dxsize, 0, dxsize, dysize );
	setLookAt (36.0,25.0,74.0, 
#ifndef __KEEP_ASPECT__
		16.0,0.0,54.0);
#else
		8.0,8.0,42.0);  // center of the house
#endif
	setFrustum (-12, 13, -2, 23, 27.75, 100);
	renderScene();

 	glutSwapBuffers();
}

void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if(h == 0)
		h = 1;

	ratio = 1.0* w / h;
	
	xsize = w;
	ysize = h;
	glViewport (0, 0, w, h);
}

int main(int argc, char **argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(xsize, ysize);
	glutCreateWindow("Perspectives With Houses");
	glutDisplayFunc(redraw);
	glutReshapeFunc(changeSize);
	// enable depth testing
	glEnable(GL_DEPTH_TEST);
	glClearColor (0.0, 0.0, 0.0, 0.0);
	glColor3f ( 1.0, 1.0, 1.0 );
	glutMainLoop();
	return 0;
}
