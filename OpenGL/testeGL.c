/*
*  testeGL.c
*  A very simple program using open GL.
*  Author : Paulo Roma Cavalcanti on 16/04/2000.
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>

typedef struct _worldPoint {
     float x, y;
} WorldPoint;

/*=========================  polyLine  ===================================*/

void polyLine ( int n, WorldPoint* pts )
{
 int i;
 glBegin ( GL_LINE_STRIP );
 for ( i = 0; i < n; ++i )
       glVertex2f ( pts[i].x, pts[i].y );
 glEnd();
}


/*=========================  polygon  ===================================*/

void polygon ( int n, WorldPoint* pts )
{
 int i;
 glBegin ( GL_LINE_LOOP );
 for ( i = 0; i < n; ++i )
       glVertex2f ( pts[i].x, pts[i].y );
 glEnd();
}


/*=========================  fillArea  ===================================*/

void fillArea ( int n, WorldPoint* pts )
{
 int i;
 glBegin ( GL_POLYGON );
 for ( i = 0; i < n; ++i )
       glVertex2f ( pts[i].x, pts[i].y );
 glEnd();
}


/*=========================  polyMarker  ===================================*/

void polyMarker ( int n, WorldPoint* pt )
{
 int i;
 glBegin ( GL_POINTS );
 for ( i = 0; i < n; ++i )
       glVertex2f ( pt[i].x, pt[i].y );
 glEnd();
}


/*=========================  setMarkerType  ===================================*/

void setMarkerType ( int type, const char* str )
{
 glPointSize ( type );
}

/*=========================  text  ===================================*/

void text ( WorldPoint* pt, const char* s )
{
 int i, size = strlen ( s );

 glRasterPos2f ( pt->x, pt->y );
 for ( i = 0; i < size; ++i )
       glutBitmapCharacter ( GLUT_BITMAP_TIMES_ROMAN_24, s[i] );

}

int xsize = 256, ysize = 256;

/*=========================  displayStuff  ===========================================*/

void displayStuff (  )
{
 WorldPoint pts1 = {0.0,0.0}, 
            pts2 = {0.0,3.0}, 
            pts3 = {4.0,3.0}, 
            pts4 = {6.0,1.5}, 
            pts5 = {4.0,0.0};
 WorldPoint pts[] = { pts1, pts2, pts3, pts4, pts5 };
 WorldPoint tori = { 0.0, -0.8 };
 WorldPoint p1 = {-1.0, -1.0}, p2 = {7.0, 4.0}; 
 int dxsize = 0.5*xsize, dysize = 0.5*ysize;
 
 setMarkerType ( 4, "+" );
 glColor3f ( 0.0, 0.0, 0.0 );
#if 0
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 gluOrtho2D ( p1.x, p2.x, p1.y, p2.y );
#else
 double fovy, z, dy, dx;
 fovy = 40.0;
 dy = (p2.y - p1.y)/2.0;
 dx = (p2.x - p1.x)/2.0;
 z = dy / tan ( fovy*3.1415927/180.0 / 2.0 );
 glMatrixMode(GL_MODELVIEW);
 glLoadIdentity();
 glTranslatef (-(p2.x-dx), -(p2.y-dy), -z);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 gluPerspective(fovy, dx/dy, z-1.0, z+1.0);
#endif
 glViewport ( 0, 0, dxsize, dysize );
 polyLine ( 5, pts );
 text ( &tori, "polyLine" );

 glViewport ( dxsize, 0, dxsize, dysize );
 polygon ( 5, pts );
 text ( &tori, "polygon" );

 glViewport ( dxsize, dysize, dxsize, dysize );
 fillArea ( 5, pts );
 text ( &tori, "fillArea" );

 glViewport ( 0, dysize, dxsize, dysize );
 polyMarker ( 5, pts );
 text ( &tori, "polyMarker" );
}


/*=========================  init  ===========================================*/

void init (  void )
{ 
 glClearColor (1.0, 1.0, 1.0, 0.0);
}


/*=========================  redraw  ===========================================*/

void redraw (  )
{ 
 glClear (GL_COLOR_BUFFER_BIT);

 displayStuff();
 
 glutSwapBuffers();
}


/*=========================  reshape  ===========================================*/

void reshape ( int w, int h )
{ 
 xsize = w;
 ysize = h;
 glViewport ( 0, 0, w, h );
}


/*=========================  keyboardHandler  ===========================================*/

void printHelp ( void )
{
 printf ( "If you went to the classes you do not need help.\n" );
}


/*=========================  keyboardHandler  ===========================================*/

void keyboardHandler ( unsigned char key, int x, int y )
{
 switch ( key )
   {
    case 'h' : printHelp(); break;
    case 'q' : exit ( 0 ); break;
   }
}

/*=========================  mouseHandler  ===========================================*/

void mouseHandler ( int button, int state, int x, int y )
{
 switch ( button )
   {
    case GLUT_LEFT_BUTTON : break;
    case GLUT_MIDDLE_BUTTON : break;
    case GLUT_RIGHT_BUTTON : break;
   }
}


/*=========================  main  ===========================================*/

int main ( int argc, char *argv[] )
{
 glutInit ( &argc, argv );
 glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGBA );
 glutInitWindowPosition ( 0, 0 );
 glutInitWindowSize ( xsize, ysize );
 glutCreateWindow ( "This is a title" );

 init();

 glutMouseFunc ( mouseHandler );
 glutKeyboardFunc ( keyboardHandler );
 glutDisplayFunc ( redraw );
 glutReshapeFunc ( reshape );
        
 glutMainLoop ( );

 return 0;
}
