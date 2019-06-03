/** @file OpenGL/drawBird.c
 *
 *  Flapping wings: CMS427, midterm exam, problem 4, page 2.
 *
 *  @author Paulo Roma Cavalcanti.
 *  @date 27/08/2017.
 *  @see http://orion.lcg.ufrj.br/cg/downloads/Computer%20Graphics%20-%20CMSC%20427.pdf
 */

#include<GL/glut.h>
#include<stdbool.h>
#include<unistd.h>

/// side of the wing.
typedef enum {LEFT=0, RIGHT=1} _side;

/// Reshape callback function.
void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (one can't make a window of zero width).
	if(h == 0)
           h = 1;

	float ratio = 1.0 * w / h;

	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45,ratio,1,1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(10,-10,10,0,0,0,0,0,1);
}

/** Draw a single triangle representing one wing of the bird.
 *
 *  @param side indicates the side of the wing for choosing its color.
 */
void drawWing(_side side) {

	if ( side == LEFT )
		glColor3f(1.0, 0, 0);
	else
		glColor3f(0, 0.5, 0.5);

	glBegin(GL_TRIANGLES);
		glVertex3f(0.0,0.0,0.0);
		glVertex3f(0.0,2.0,0.0);
		glVertex3f(3.0,0.0,0.0);
	glEnd();

}

/** Draw a single wing of the bird.
 *
 *  Since the fixed point of all transformations is the origin, 
 *  it is easier to rotate first about the common axis (x axis at this moment), 
 *  by the angle controlling the animation. 
 *
 *  Then, rotate by -90 degrees (clockwise direction) and scale
 *  by a factor of 1.5, or -1.5 for the left wing (to reflect it about the y axis). 
 *  The last step is to translate by (2,3) to the final position.
 *
 *  @param side indicates the side of the wing: left or right.
 *  @param angle between wings.
 */
void drawBird(_side side, float angle) {

	glPushMatrix ();

	glTranslatef (2.0, 3.0, 0.0);
	glScalef (side==LEFT?-1.5:1.5, 1.5, 1.5);
	glRotatef (-90.0, 0.0, 0.0, 1.0);
	glRotatef (angle, 1.0, 0.0, 0.0);

	drawWing(side);

	glPopMatrix();
}

/** Render the whole scene, by drawing the left and right wings.
 *
 *  The rotation angle should be incremented/decremented 
 *  (to create the sensation of flapping wings)
 *  accordingly, to keep it in the range [0,60] degrees.
 */
void renderScene(void) {
	// rotation angle.
	static float angle = 0.0;

	// controls the direction of the animation, up or down.
	static bool down = false;

	// clear the window, and the depth buffer.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	drawBird(LEFT, angle);
	drawBird(RIGHT, angle);

	if (angle >= 60) {
		down = true;
	}
	else if (angle <= 0) {
		down = false;
	}
	down ?  angle-- : angle++;

	glutSwapBuffers();
}

int main(int argc, char **argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(420,420);
	glutCreateWindow("Flapping wings: CMS427, midterm exam, problem 4, page 2");
	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutMainLoop();
	return 0;
}

