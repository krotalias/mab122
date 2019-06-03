/**
 *  @file torus.cpp
 *
 *  This program demonstrates the creation of display lists for constructing some quadric surfaces.
 *
 *  Quadric surfaces are defined by quadratic equations in two dimensional space.
 *
 *  Display list is a group of OpenGL commands that have been stored (compiled) for later execution. 
 *  Once a display list is created, all vertex and pixel data are evaluated and copied into the 
 *  display list memory on the server machine. It is only one time process. After the display list has been prepared (compiled), 
 *  you can reuse it repeatedly without re-evaluating and re-transmitting data over and over again to draw each frame. 
 *  Display list is one of the fastest methods to draw static data because vertex data and OpenGL commands 
 *  are cached in the display list and minimize data transmissions from the client to the server side. 
 *  It means that it reduces CPU cycles to perform the actual data transfer.
 *
 *  - To Compile:
 *    - gcc torus.cpp -o torus-linux64
 *  - or, for MacOS:   
 *    - gcc torus.cpp -o torus.osx
 *
 *  - Usage: 
 *  <PRE>
 *     - torus.osx R   r   r'  nR  nr
 *     - torus.osx 1.0 0.3 0.5 25  10
 *  </PRE>
 *   
 *  - Requirements:
 *    - DevIL-ILUT-devel (Fedora)
 *    - libdevil1c2 libdevil-dev (Ubuntu)
 *    - libdevil (MacPorts)
 *
 *  @author Paulo R Cavalcanti
 *  @date 06/11/2017
 *  @see https://www.opengl.org/archives/resources/code/samples/redbook/torus.c
 *  @see http://www.songho.ca/opengl/gl_displaylist.html
 *  @see https://web.cs.wpi.edu/~matt/courses/cs563/talks/renderman/quadric.html
 *  @see http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4037402
 *  @see http://www.nptel.ac.in/courses/Webcourse-contents/IIT-Delhi/Computer%20Aided%20Design%20&%20ManufacturingI/mod3/01.htm
 *  @see http://openil.sourceforge.net/
 *  @see http://tutorial.math.lamar.edu/Classes/CalcIII/QuadricSurfaces.aspx
 *  @see http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/
 */

#define GL_GLEXT_PROTOTYPES

#include <IL/il.h>
#include <GL/glut.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <Arcball/Arcball.hpp>

/// &Pi; definition. Unused. We are using M_PI from math.h
/// The difference between const and define
/// is that const has a type and define is trick,
/// because it is the pre-compiler which replaces the value.
#define M_PI_ 3.14159265358979323846
/// Texture Image width.
#define checkImageWidth 128
/// Texture Image height.
#define checkImageHeight 128
/// Far plane distance.
/// Near and far clipping planes should always be specified as positive. 
/// That's because they represent how far in front of the "camera" the clipping planes will be.
#define farPlane (zobs+10*d)
/// Front plane distance.
#define frontPlane 1
/// Return a positive angle.
#define posAng(x) ((x) < 0.0 ? ((x)+TWOPI) : (x))
/// Two PI.
#define TWOPI (2*M_PI)
/// Convert an angle to radians.
#define toRad(x) ((x)*M_PI/180.0)
/// Convert an angle to degrees.
#define toDeg(x) ((x)*180.0/M_PI)
/// Clamp to [0,1].
#define clamp(x) (fmin(fmax((x),0),1))
/// Vector length.
#define LEN(x,y,z) (sqrt((x)*(x)+(y)*(y)+(z)*(z)))
/// Axes.
#define theAxes (2*nObjects)
/// Number of objects;
#define nObjects (sizeof(objNames)/sizeof(char*))
/// Scale for drawing normals.
#define scale (0.05*diag)  

/// Camera position. 
double zobs = 10.0;
/// Camera viewing angle, also known as opening angle.
double viewingAngle = 45;
/// Object display list base index.
GLuint objects;
/// Distance from the center of the tube to the center of the torus.
double R1 = 1.0;
/// Radius of the tube.
double R2 = 0.3;
/// Radius in z direction for an ellipsoid.
double R3 = 0.5;
/// Number of samples.
int nR1 = 25;
/// Number of samples.
int nR2 = 10;
/// Whether draw wireframe or filled polygons.
int drawWire = 0;
/// Whether display vertex normals.
int drawNormals = 0;
/// Whether to apply texture.
int applyTexture = 0;
/// Whether to show the coordinate axes.
int showAxes = 0;
/// Whether to show the bounding boxes.
int showBox = 0;
/// Whether to use a procedural texture.
int procImage = 1;
/// Whether to use the arcball paradigm.
int useArcball = 0;
/// Automatic spin.
int spin = 1;
/// 3D point structure.
typedef struct _Point {
  double x, y, z;
} Point;
/// Axis aligned Bounding box.
typedef Point _BBOX[2];

/// Quadric types. 
typedef enum _oType {_torus = 0, _cylinder, _cone, _ellipsoid, _paraboloid, _hyperbolic_paraboloid, _hyperboloid} oType;

/// Rotation axes.
typedef enum _rAxes {xAxis, yAxis, zAxis} rAxes;

/// Current axis.
rAxes axis = xAxis;

/// Rotation angle increment.
GLfloat rotAngle = 2.0;

/// Controls the trigger of the timer function.
unsigned int delay = 60;

/// Procedural image for a checkerboard.
GLubyte checkImage[checkImageHeight][checkImageWidth][4];
/// Texture identifier.
GLuint texName;
/// IL image identifier.
ILuint ilImage;
/// Object names.
const char *objNames[] = {"Torus", "Cylinder", "Cone", "Ellipsoid", "Paraboloid", "Hyperbolic Paraboloid", "Hyperboloid"};
/// Object boxes.
_BBOX objBoxes[nObjects];
/// Current object.
int OBJECT = _torus;
/// ArcBall.
Arcball arcball;

/** Add a new point to a bounding box.
 *
 * The axis-aligned minimum bounding box (or AABB) for a given point set is its minimum bounding box
 * subject to the constraint that the edges of the box are parallel to the (Cartesian) coordinate axes.
 * It is simply the Cartesian product of N intervals each of which is defined by the minimal and maximal
 * value of the corresponding coordinate for the points in S.
 *
 * Axis-aligned minimal bounding boxes are used to an approximate location of an object in question and 
 * as a very simple descriptor of its shape. 
 * For example, in computational geometry and its applications when it is required to find intersections in the set of objects, 
 * the initial check is the intersections between their MBBs. Since it is usually a much less expensive operation than 
 * the check of the actual intersection (because it only requires comparisons of coordinates),
 * it allows to quickly exclude from checks the pairs that are far apart.
 *
 *  @param index Position into objBoxes array.
 *  @param x Point coordinate.
 *  @param y Point coordinate.
 *  @param z Point coordinate.
 *  @see https://en.wikipedia.org/wiki/Minimum_bounding_box
 */
void updateBox(oType index, double x, double y, double z) {
   _BBOX *bbox = objBoxes+index;

   (*bbox)[0].x = fmin(x,(*bbox)[0].x);
   (*bbox)[1].x = fmax(x,(*bbox)[1].x);
   (*bbox)[0].y = fmin(y,(*bbox)[0].y);
   (*bbox)[1].y = fmax(y,(*bbox)[1].y);
   (*bbox)[0].z = fmin(z,(*bbox)[0].z);
   (*bbox)[1].z = fmax(z,(*bbox)[1].z);
}

/** Return bounding box dimension.
 *
 *  @param index Position into objBoxes array.
 *  @param dx Bounding box length.
 *  @param dy Bounding box width.
 *  @param dz Bounding box height.
 *  @return Bounding box diameter.
 */
double getBoxSize(oType index, double *dx, double *dy, double *dz) {
   _BBOX *bbox = objBoxes+index;

   *dx = (*bbox)[1].x - (*bbox)[0].x;
   *dy = (*bbox)[1].y - (*bbox)[0].y;
   *dz = (*bbox)[1].z - (*bbox)[0].z;

   return LEN(*dx,*dy,*dz);
}

/** Draw a cone, centered at the origin, with a circular base on plane XY.
 *
 * Implicit Equation:
 * - @f${x^{2} \over c^{2}}+{y^{2} \over c^{2}} = {z^{2}},\ c = \frac{r}{h}@f$
 * - 
 * - @f${f(x,y,z) = \frac{h}{r} x^{2}}+{\frac{h}{r} y^{2}}-{\frac{r}{h} z^{2} = 0}@f$
 * - 
 * - @f$\nabla f(x,y,z) = (\frac{h}{r} x,\ \frac{h}{r} y,\ -{\frac{r}{h} z})@f$
 *
 * Parametric Equation:
 * - x(u,v) = u cos(v), 
 * - y(u,v) = u sin(v), 
 * - z(u,v) = u * h / r, 
 * - u &isin; [0,h], v &isin; [0,2&Pi;]
 * - df/du = (cos(v), sin(v), h/r), 
 * - df/dv = (-u sin(v), u cos(v) , 0)
 * - df/du x df/dv = (h/r * u cos(v), h/r * u sin(v), -u) = (h/r * x(u,v), h/r * y(u,v), -r/h * z(u,v)) 
 *
 * Inverse parametrization:
 * - u = (r/h) * z, or u = (r/h) * (h-z), if cone is upside down.
 * - v = atan2(y,x), if v < 0 then v += 2&Pi;
 *
 * Note: atan2 returns values &isin; [-&Pi;, &Pi;]
 *
 * @param h cone height.
 * @param r cone base radius.
 * @param vs number of divisions in z direction.
 * @param rs number of angular divisions.
 * @param normals whether to draw normals or polygons.
 * @see http://mathinsight.org/parametrized_surface_examples
 * @see http://mathworld.wolfram.com/Cone.html
 * @see https://en.wikipedia.org/wiki/Atan2
 */
void cone(double h, double r, int vs, int rs, int normals) {
   int i, j, k;
   double z, th, x, y, ln, nx, ny, nz, u;
   double hr = h / r;

   double diag;
   if (normals) {
      double dx,dy,dz;
      diag = getBoxSize(_cone,&dx,&dy,&dz);
   }

   for (i = 0; i < vs; i++) {
      glBegin(normals?GL_LINES:GL_QUAD_STRIP);
      for (j = 0; j <= rs; j++) {
         for (k = 1; k >= 0; k--) {
            z = h * (i+k) / vs;              
            u = (h-z) / hr; // cone is upside down

            th = j * TWOPI / rs;
            x = u * cos(th);
            y = u * sin(th);

            nx = x * hr; ny = y * hr; nz = u; // cone is upside down
            ln = LEN(nx,ny,nz);
            nx /= ln; ny /= ln; nz /= ln;
            th = clamp(th/TWOPI);

            glNormal3d(nx, ny, nz); 
            glTexCoord2d(th,z/h); 
            glVertex3d(x, y, z);
            if (normals) 
               glVertex3d(x+scale*nx, y+scale*ny, z+scale*nz);
            else
               updateBox(_cone, x, y, z);
         }
      }
      glEnd();
   }
}

/** Draw a cylinder, centered at the origin, with a circular base on plane XY.
 *
 * Implicit Equation:
 * - @f${x^{2}+y^{2}} = {r^{2}},@f$ any z
 * - 
 * - @f${f(x,y,z) = x^{2}+y^{2}}-{r^{2} = 0}@f$
 * -
 * - @f${\nabla f(x,y,z) = (x, y, 0)}@f$
 *
 * Parametric Equation:
 * - x(u,v) = r cos(v), 
 * - y(u,v) = r sin(v), 
 * - z(u,v) = u, 
 * - u &isin; [0,h], v &isin; [0,2&Pi;]
 * - df/du = (0, 0, 1), 
 * - df/dv = (-r sin(v), r cos(v) , 0)
 * - df/dv x df/du = (r cos(v), r sin(v), 0) = (x(u,v), y(u,v), 0) 
 *
 * Inverse parametrization:
 * - u = z
 * - v = acos(x/r) = atan2(y,x), if v < 0 then v += 2&Pi;
 *
 * @param h cylinder height.
 * @param r cylinder radius.
 * @param vs number of divisions in z direction.
 * @param rs number of angular divisions.
 * @param normals whether to draw normals or polygons.
 * @see http://mathworld.wolfram.com/Cylinder.html
 * @see https://en.wikipedia.org/wiki/Quadric
 * @see http://mathworld.wolfram.com/QuadraticSurface.html
 */
void cylinder(double h, double r, int vs, int rs, int normals) {
   int i, j;
   double z0, z1, th, x, y, ln, nx, ny;

   double diag;

   if (normals) {
      double dx,dy,dz;
      diag = getBoxSize(_cylinder,&dx,&dy,&dz);
   }

   for (i = 0; i < vs; i++) {       // generate stacks
      z0 = h * i / vs;              // stack bottom z
      z1 = h * (i+1) / vs;          // stack top z
      glBegin(normals?GL_LINES:GL_QUAD_STRIP);
      for (j = 0; j <= rs; j++) {   // generate slices
           th = j * TWOPI / rs;     // slice left theta
           x = r * cos(th);
           y = r * sin(th);

           ln = LEN(x,y,0);
           nx = x/ln; ny = y/ln;
           th = clamp(th/TWOPI);

           glNormal3d(nx, ny, 0); 
           glTexCoord2d(th,z1/h); 
           glVertex3d(x, y, z1);
           if (normals) 
               glVertex3d(x+scale*nx, y+scale*ny, z1);
           else
               updateBox(_cylinder, x, y, z1);
           glNormal3d(nx, ny, 0); 
           glTexCoord2d(th,z0/h); 
           glVertex3d(x, y, z0);
           if (normals) 
              glVertex3d(x+scale*nx, y+scale*ny, z0);
           else
              updateBox(_cylinder, x, y, z0);
      }
      glEnd();
   }
}

/** Draw a paraboloid, centered at the origin, with a circular base on plane XY.
 *
 * Implicit Equation:
 * - @f${x^{2}+y^{2}} = {c\ z},\ c = \frac{r^2}{h}@f$
 * -
 * - @f${f(x,y,z) = \frac{h}{r^2} x^{2} + \frac{h}{r^2} y^{2} - z = 0}@f$
 * -
 * - @f${\nabla f(x,y,z) = (2 \frac{h}{r^2} x,\ 2 \frac{h}{r^2} y,\ - 1)}@f$
 *
 * Parametric Equation:
 * - x(u,v) = u cos(v), 
 * - y(u,v) = u sin(v), 
 * - z(u,v) = h * (u/r)², (z = h -> u = r)
 * - u &isin; [0,h], v &isin; [0,2&Pi;]
 * - df/du = (cos(v), sin(v), 2*u h/r²), 
 * - df/dv = (-u sin(v), u cos(v) , 0)
 * - df/dv x df/du = (2*u h/r² u cos(v), 2*u h/r² u sin(v), -u) = (2h/r² x(u,v), 2h/r² y(u,v), -1) 
 *
 * Inverse parametrization:
 * - u = sqrt(z * r² / h ) or u = sqrt( (h-z) * r² / h ), if paraboloid is upside down.
 * - v = atan2(y,x), if v < 0 then v += 2&Pi;
 *
 * @param h paraboloid height.
 * @param r paraboloid radius.
 * @param vs number of divisions in z direction.
 * @param rs number of angular divisions.
 * @param normals whether to draw normals or polygons.
 * @see http://mathworld.wolfram.com/Paraboloid.html
 * @see https://en.wikipedia.org/wiki/Quadric
 * @see http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node61.html
 */
void paraboloid(double h, double r, int vs, int rs, int normals) {
   int i, j, k;
   double z, th, x, y, ln, nx, ny, nz, u;
   double hr2 = h/(r*r);

   double diag;
   if (normals) {
      double dx,dy,dz;
      diag = getBoxSize(_paraboloid,&dx,&dy,&dz);
   }

   for (i = 0; i < vs; i++) {
      glBegin(normals?GL_LINES:GL_QUAD_STRIP);
      for (j = 0; j <= rs; j++) {
         for (k = 1; k >= 0; k--) {
            z = h * (i+k) / vs;              
            u = sqrt((h-z)/hr2); // paraboloid is upside down

            th = j * TWOPI / rs;
            x = u * cos(th);
            y = u * sin(th);

            nx = 2 * x * hr2; ny = 2 * y * hr2; nz = 1; // paraboloid is upside down
            ln = LEN(nx,ny,nz);
            nx /= ln; ny /= ln; nz /= ln;
            th = clamp(th/TWOPI);

            glNormal3d(nx, ny, nz); 
            glTexCoord2d(th,z/h); 
            glVertex3d(x, y, z);
            if (normals) 
               glVertex3d(x+scale*nx, y+scale*ny, z+scale*nz);
            else
               updateBox(_paraboloid, x, y, z);
         }
      }
      glEnd();
   }
}

/** Draw a one-sheeted hyperboloid of revolution, centered at the origin, with a circular base parallel to plane XY.
 *
 * Implicit Equation:
 * - @f${{x^{2} \over a^{2}} + {y^{2} \over a^{2}} - {z^{2} \over c^{2}}} = 1@f$
 * -
 * - @f$r = a * \sqrt{1 + \frac{h^2}{4c^2}}@f$
 * -  
 * - @f$c = \sqrt{h^2 / (4 (\frac{r^2}{a^2} - 1))},\ r \ge a@f$
 * -
 * - @f${f(x,y,z) = {c^2} x^{2} + {c^2} y^{2} - {a^2} z^{2} - {a^2} {c^2} = 0}@f$
 * -
 * - @f${\nabla f(x,y,z) = (2 \frac{x}{a^2},\ 2 \frac{y}{a^2},\ -2 \frac{z}{c^2})}@f$
 *
 * Parametric Equation:
 * - x(u,v) = a sqrt(1+u²) cos(v), 
 * - y(u,v) = a sqrt(1+u²) sin(v), 
 * - z(u,v) = c u
 * - u &isin; [-h,h], v &isin; [0,2&Pi;]
 * - df/du = (-a u/srqt(1+u²) cos(v), -a u/srqt(1+u²) sin(v), c), 
 * - df/dv = (-a sqrt(1+u²) sin(v), a sqrt(1+u²) cos(v) , 0)
 * - df/dv x df/du = (c x, c y, -a²/c z)
 *
 * Inverse parametrization:
 * - u = sqrt( x² / a² + y² / a² - 1 )
 * - v = atan2(y,x), if v < 0 then v += 2&Pi;
 *
 * @param h hyperboloid height.
 * @param a hyperboloid skirt radius.
 * @param r hyperboloid base radius.
 * @param vs number of divisions in z direction.
 * @param rs number of angular divisions.
 * @param normals whether to draw normals or polygons.
 * @see http://mathworld.wolfram.com/One-SheetedHyperboloid.html
 * @see http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node61.html
 * @see https://rechneronline.de/pi/hyperboloid-e.php
 */
void hyperboloid(double h, double a, double r, int vs, int rs, int normals) {
   int i, j, k;
   double z, th, x, y, ln, nx, ny, nz, u;
   double c = sqrt( h*h / (4 * (r*r/(a*a) - 1)) );
   if ( r < a ) {
      printf ( "Hyperboloid: r < a -> invalid!\n" );
      return;
   }

   double diag;
   if (normals) {
      double dx,dy,dz;
      diag = getBoxSize(_hyperboloid,&dx,&dy,&dz);
   }
   
   for (i = -vs; i < vs; i++) {
      glBegin(normals?GL_LINES:GL_QUAD_STRIP);
      for (j = 0; j <= rs; j++) {
         for (k = 1; k >= 0; k--) {
            z = h * (i+k) / vs;              
            u = z/c;
            u = a * sqrt(1+u*u);

            th = j * TWOPI / rs;
            x = u * cos(th);
            y = u * sin(th);

            nx = x * c; ny = y * c; nz = -z * a * a / c;
            ln = LEN(nx,ny,nz);
            nx /= ln; ny /= ln; nz /= ln;
            th = clamp(th/TWOPI);

            glNormal3d(nx, ny, nz); 
            glTexCoord2d(th,(z+h)/(2*h)); // map [0,1] to [-h,h]
            glVertex3d(x, y, z);
            if (normals) 
               glVertex3d(x+scale*nx, y+scale*ny, z+scale*nz);
            else
               updateBox(_hyperboloid, x, y, z);
         }
      }
      glEnd();
   }
}

/** Draw a hyperbolic paraboloid.
 *
 * Implicit Equation:
 * - @f${xy = z}@f$
 * -
 * - @f${f(x,y,z) = xy - z = 0}@f$
 * -
 * - @f${\nabla f(x,y,z) = (y, x, -1)}@f$
 *
 * Parametric Equation:
 * - x(u,v) = u, 
 * - y(u,v) = v, 
 * - z(u,v) = uv 
 * - u &isin; [0,w], v &isin; [0,h]
 * - df/du = (1, 0, v), 
 * - df/dv = (0, 1, u)
 * - df/du x df/dv = (y, x, -1) 
 *
 * Inverse parametrization:
 * - u = x
 * - v = y
 *
 * @param w hyperbolic paraboloid width.
 * @param h hyperbolic paraboloid length.
 * @param ws number of divisions in x direction.
 * @param hs number of divisions in y direction.
 * @param normals whether to draw normals or polygons.
 * @see http://mathworld.wolfram.com/HyperbolicParaboloid.html
 * @see https://en.wikipedia.org/wiki/Quadric
 * @see http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node61.html
 */
void hyperbolic_paraboloid(double w, double h, int ws, int hs, int normals) {
   int i, j, k;
   double u, v,  z, x, y, ln, nx, ny, nz;

   double diag;
   if (normals) {
      double dx,dy,dz;
      diag = getBoxSize(_hyperbolic_paraboloid,&dx,&dy,&dz);
   }
  
   for (i = -hs; i < hs; i++) {
      glBegin(normals?GL_LINES:GL_QUAD_STRIP);
      for (j = -ws; j <= ws; j++) {
         for (k = 1; k >= 0; k--) {
            u = j*w/ws;
            v = (i+k)*h/hs;	 

            x = u;
            y = v;
            z = u*v;

            nx = y; ny = x; nz = -1; 
            ln = LEN(nx,ny,nz);
            nx /= ln; ny /= ln; nz /= ln;

            glNormal3d(nx, ny, nz); 
            glTexCoord2d((u+w)/(2*w),(v+h)/(2*h)); // map [0,1] to [-h,h]
            glVertex3d(x, y, z);
            if (normals) 
               glVertex3d(x+scale*nx, y+scale*ny, z+scale*nz);
            else
               updateBox(_hyperbolic_paraboloid, x, y, z);
         }
      }
      glEnd();
   }
}

/** Draw a torus, centered at the origin, and azimuthally symmetric about the z-axis.
 *
 * In geometry, a torus is a surface of revolution generated by revolving a circle 
 * in three-dimensional space about an axis coplanar with the circle.
 * - Volume: 2 × π² × R × r²
 * - Surface area: 4 × π² × R × r
 * - Common objects with this shape: Doughnut, Ring, Lifebuoy
 *
 * Implicit Equation:
 * - @f$(r_1-\sqrt{(x^2+y^2)})^2+z^2=r_2^2,@f$
 * -
 * - @f$f(x,y,z) = (x^2+y^2+z^2+r_1^2-r_2^2)^2 - 4 r_1^2 (x^2+y^2) = 0@f$, a quartic equation
 *
 * Parametric Equation:
 * - x(u,v) = (r1 + r2 cos(v)) cos(u), 
 * - y(u,v) = (r1 + r2 cos(v)) sin(u), 
 * - z(u,v) = r2 sin(v), 
 * - u &isin; [0,2&Pi;], v &isin; [0,2&Pi;]
 * - df/du = (-y(u,v), x(u,v), 0), 
 * - df/dv = (-z(u,v) cos(u), -z(u,v) sin(u), r2 cos(v))
 * - df/du x df/dv = (r2 x(u,v) cos(v), r2 y(u,v) cos(v), z(u,v) (y(u,v) sin(u) + x(u,v) (cos(u))) = (r2 cos(v) x(u,v), r2 cos(v) y(u,v), (r1 + r2 cos(v)) z(u,v))
 * - df/du x df/dv = (cos(v) cos(u), cos(v) sin(u), sin(v))
 *
 * Inverse parametrization:
 * - u = atan2(y,x), if u < 0 then u += 2&Pi;
 * - v = asin(z/r2)
 *
 * When:
 * - r1 > r2, the surface will be the familiar ring torus.
 * - r1 = r2 corresponds to the horn torus, which in effect is a torus with no "hole".
 * - r1 < r2 describes the self-intersecting spindle torus.
 * - r1 = 0, the torus degenerates to the sphere.
 *
 * @param r1 distance from the center of the tube to the center of the torus.
 * @param r2 radius of the tube.
 * @param numc number of divisions of the circle with radius r2.
 * @param numt number of divisions of the circle with radius r1.
 * @param normals whether to draw normals or polygons.
 * @see http://mathworld.wolfram.com/Torus.html
 * @see https://en.wikipedia.org/wiki/Torus
 * @see https://blogs.scientificamerican.com/roots-of-unity/a-few-of-my-favorite-spaces-the-torus/
 * @see https://www.mathsisfun.com/geometry/torus.html
 * @see http://trecs.se/torus.php
 * @see http://mathworld.wolfram.com/QuarticEquation.html
 */
void torus(double r1, double r2, int numc, int numt, int normals) {
   int i, j, k;
   double s, t, x, y, z;
   double u, v, cu, cv, su, sv, nx, ny, nz;
   double ln;

   double diag;
   if (normals) {
      double dx,dy,dz;
      diag = getBoxSize(_torus,&dx,&dy,&dz);
   }

   for (i = 0; i < numc; i++) {
      glBegin(normals?GL_LINES:GL_QUAD_STRIP);
      for (j = 0; j <= numt; j++) {
         for (k = 1; k >= 0; k--) {
#if 0           
            s = (i + k) % numc + 0.5;
            t = j % numt;
#else
            s = (i + k); 
            t = j;
#endif
            u = t*TWOPI/numt;
            v = s*TWOPI/numc;

            cu = cos(u);
            cv = cos(v);
            su = sin(u);
            sv = sin(v);

            // (x,y,z) = f(u,v)
            x = (r1 + r2 * cv) * cu;
            y = (r1 + r2 * cv) * su;
            z = r2 * sv;

            // df/du x df/dv
            nx = cv * cu;
            ny = cv * su;
            nz = sv;

            ln = LEN(nx,ny,nz);
            nx /= ln; ny /= ln; nz /= ln;

            u = clamp(u/TWOPI);
            v = clamp(v/TWOPI);

            glTexCoord2d(u,v);

            glNormal3d(nx, ny, nz);
            glVertex3d(x, y, z);
            if (normals) 
               glVertex3d(x+scale*nx, y+scale*ny, z+scale*nz);
            else
               updateBox(_torus, x, y, z);
         }
      }
      glEnd();
   }
}

/** Draw an ellipsoid or a sphere centered at the origin.
 *
 * Implicit Equation:
 * - @f${x^{2} \over a^{2}}+{y^{2} \over b^{2}}+{z^{2} \over c^{2}}=1,\ a= r_1,\ b = r_2,\ c = r_3@f$
 * -
 * - @f${f(x,y,z) = r_2^2 r_3^2\ x^{2}+r_1^2 r_3^2\ y^{2}+r_1^2 r_2^2\ z^{2}} - (r_1 r_2 r_3)^2 = 0@f$
 * -
 * - @f${\nabla f(x,y,z) = (r_2^2 r_3^2\ x,\ r_1^2  r_3^2\ y,\ r_1^2  r_2^2\ z)}@f$
 *
 * A sphere is defined as the set of all points in three-dimensional Euclidean space R³
 * that are located at a distance r (the "radius") from a given point (the "center"). 
 *
 * An ellipsoid is a surface that may be obtained from a sphere by deforming it 
 * by means of directional scalings, or more generally, of an affine transformation.
 *
 * Parametric Equation:
 * - x(u,v) = r1 cos(u) sin(v), 
 * - y(u,v) = r2 sin(u) sin(v), 
 * - z(u,v) = r3 cos(v), 
 * - u &isin; [0,2&Pi;], v &isin; [0,&Pi;]
 * - df/du = (-r1 sin(u) sin(v), r2 cos(u) sin(v), 0), 
 * - df/dv = (r1 cos(u) cos(v), r2 sin(u) cos (v), -r3 sin(v))
 * - df/dv x df/du = (r2 r3/r1 x(u,v), r1 r3/r2 y(u,v), r1 r2/r3 z(u,v)) sin v = (r2² r3² x(u,v), r1² r3² y(u,v), r1² r2² z(u,v))
 *
 * Inverse parametrization:
 * - u = atan2((y,x) * (r1/r2)), if u < 0 then u += 2&Pi;
 * - v = acos(z/r3)
 *
 * @param r1 x radius of the ellipsoid.
 * @param r2 y radius of the ellipsoid.
 * @param r3 z radius of the ellipsoid.
 * @param numc number of divisions along the zero degree longitude (number of latitude divisions).
 * @param numt number of divisions along the zero degree latitude - equator (number of longitude divisions).
 * @param normals whether to draw normals or polygons.
 * @see http://mathworld.wolfram.com/Ellipsoid.html
 * @see http://mathworld.wolfram.com/Sphere.html
 * @see https://en.wikipedia.org/wiki/Ellipsoid
 * @see https://sciencing.com/equators-latitude-6314100.html
 */
void sphere(double r1, double r2, double r3, int numc, int numt, int normals) {
   int i, j, k;
   double s, x, y, z;
   double u, v, cu, cv, su, sv, nx, ny, nz;
   double ln;

   double diag;
   if (normals) {
      double dx,dy,dz;
      diag = getBoxSize(_ellipsoid,&dx,&dy,&dz);
   }

   for (i = 0; i < numc; i++) {
      glBegin(normals?GL_LINES:GL_QUAD_STRIP);
      for (j = 0; j <= numt; j++) {
         for (k = 1; k >= 0; k--) {
            s = (i + k); 

            u = j*TWOPI/numt;
            v = s*M_PI/numc;

            cu = cos(u);
            cv = cos(v);
            su = sin(u);
            sv = sin(v);

            // (x,y,z) = f(u,v)
            x = r1 * cu * sv;
            y = r2 * su * sv;
            z = r3 * cv;

            // df/du x df/dv
            nx = r2*r3*r2*r3 * x;
            ny = r1*r3*r1*r3 * y;
            nz = r1*r2*r1*r2 * z;

            ln = LEN(nx,ny,nz);
            nx /= ln; ny /= ln; nz /= ln;

            u = clamp(u/TWOPI);
            v = clamp(v/M_PI);

            glTexCoord2d(u,v);

            glNormal3d(nx, ny, nz);
            glVertex3d(x, y, z);
            if (normals) 
               glVertex3d(x+scale*nx, y+scale*ny, z+scale*nz);
            else
               updateBox(_ellipsoid, x, y, z);
         }
      }
      glEnd();
   }
}

/** Draw the three coordinate axes.
 *
 *  The length of the axes is r1 + r2.
 *
 *  @param r1 axis length.
 *  @param r2 axis additional length.
 */
void drawAxes (double r1, double r2) {
   double d;
   d = r1+r2;
   glBegin(GL_LINES);
      glColor3d(1.0,0.0,0.0);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(d,0.0,0.0);
      glColor3d(0.0,1.0,0.0);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(0.0,d,0.0);
      glColor3d(0.0,0.0,1.0);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(0.0,0.0,d);
   glEnd();
}

/** Draw the bounding box as a dashed hexahedron.
 *
 *  @param BBOX given bounding box.
 */
void drawBox (const _BBOX BBOX) {
    // glPushAttrib is done to return everything to normal after drawing
    glPushAttrib(GL_ENABLE_BIT);
    glLineStipple(3, 0xAAAA);  // [1]
    glEnable(GL_LINE_STIPPLE);

    glBegin(GL_LINE_LOOP);
    glColor3d(1.0,1.0,1.0);
    glVertex3d(BBOX[0].x, BBOX[0].y, BBOX[0].z);
    glVertex3d(BBOX[1].x, BBOX[0].y, BBOX[0].z);
    glVertex3d(BBOX[1].x, BBOX[1].y, BBOX[0].z);
    glVertex3d(BBOX[0].x, BBOX[1].y, BBOX[0].z);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex3d(BBOX[0].x, BBOX[0].y, BBOX[1].z);
    glVertex3d(BBOX[1].x, BBOX[0].y, BBOX[1].z);
    glVertex3d(BBOX[1].x, BBOX[1].y, BBOX[1].z);
    glVertex3d(BBOX[0].x, BBOX[1].y, BBOX[1].z);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex3d(BBOX[0].x, BBOX[0].y, BBOX[0].z);
    glVertex3d(BBOX[1].x, BBOX[0].y, BBOX[0].z);
    glVertex3d(BBOX[1].x, BBOX[0].y, BBOX[1].z);
    glVertex3d(BBOX[0].x, BBOX[0].y, BBOX[1].z);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex3d(BBOX[0].x, BBOX[1].y, BBOX[0].z);
    glVertex3d(BBOX[1].x, BBOX[1].y, BBOX[0].z);
    glVertex3d(BBOX[1].x, BBOX[1].y, BBOX[1].z);
    glVertex3d(BBOX[0].x, BBOX[1].y, BBOX[1].z);
    glEnd();

    glPopAttrib();
}

/** Create checkerboard texture.
 *
 *  - If i's 4th bit is set, i.e., every 8 iterations of i, i & 0x8 yields nonzero. 
 *  - The same for j. If either is nonzero, but the other is zero, then XOR will yield 1. 
 *  - This is multiplied by 255, i.e., the maximum value of a channel.
 *  <PRE> 
 *     -  i    -> i & 1000
 *     - 0-7   -> 0
 *     - 8-15  -> 8
 *     - 16-23 -> 0
 *     - 24-31 -> 8
 *     - ....
 *  </PRE>
 *
 *  @see https://www.tutorialspoint.com/cprogramming/c_bitwise_operators.htm
 *  @see https://www.miniwebtool.com/bitwise-calculator/
 */
void makeCheckImage(void) {
   int i, j, c;

   for (i = 0; i < checkImageHeight; i++) {
      for (j = 0; j < checkImageWidth; j++) {
         // (i and 1000) xor (j and 1000)
         c = ( (i&0x8)==0 ) ^ ( (j&0x8)==0 );
         c *= 255;
         if ( (i&0x8)==0 ) {
            checkImage[i][j][0] = (GLubyte) fmax(c,255);  // red
            checkImage[i][j][1] = (GLubyte) c;            // green
         } else {
            checkImage[i][j][0] = (GLubyte) c;            // red
            checkImage[i][j][1] = (GLubyte) fmax(c,128);  // green
         }
         checkImage[i][j][2] = (GLubyte) c;               // blue
         checkImage[i][j][3] = (GLubyte) 255;             // alpha
      }
   }
}

/** Load an image from a file.
 *
 *  @param filename file name.
 *  @return -1 on fail, or image name on success.
 *  @see http://openil.sourceforge.net/tuts/tut_4/index.htm
 */
int loadImage (const char *filename) {
    ILuint    imageName;

    ilGenImages(1, &imageName);    // load just one image
    ilBindImage(imageName);        // set image name as the current image in DevIL

    // load image filename using DevIL
    if ( !ilLoadImage(filename) )
         return -1;

    // convert every colour component to unsigned byte
    // one can replace IL_RGB with IL_RGBA if the image contains alpha channel
    if ( !ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE) )
         return -1;

    return imageName;
}

/** Initialize texture stuff.
 *
 *  Create a single texture called texName.
 *
 *  @see https://gregs-blog.com/2008/01/17/opengl-texture-filter-parameters-explained/
 *  @see https://open.gl/textures
 */
void initTexture(void) {
   // Create a checkerboard, by using a procedural texture: red, green and white.
   makeCheckImage();

   // UNSIGNED_BYTE requires one byte alignment.
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   // Create a single texture and save its name as texName.
   glGenTextures(1, &texName);
   // Set texName as the current texture.
   glBindTexture(GL_TEXTURE_2D, texName);
   // Specifies a texture environment. 
   // GL_DECAL ignores the primary color and just shows the texel.
   glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
   // GL_TEXTURE_MIN_FILTER is used whenever a surface is rendered with 
   // smaller dimensions than its corresponding texture bitmap (far away objects). 
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
   // GL_TEXTURE_MAG_FILTER is used in the exact opposite case: 
   // a surface is bigger than the texture being applied (near objects)
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

   if ( procImage ) {
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, checkImageWidth, checkImageHeight, 
                0, GL_RGBA, GL_UNSIGNED_BYTE, checkImage);
   } else {

      glTexImage2D(GL_TEXTURE_2D, 0, ilGetInteger(IL_IMAGE_BPP), ilGetInteger(IL_IMAGE_WIDTH), ilGetInteger(IL_IMAGE_HEIGHT),
                0, ilGetInteger(IL_IMAGE_FORMAT), GL_UNSIGNED_BYTE, ilGetData());
   }

   // Mipmap only makes sense for the MIN_FILTER.
   glGenerateMipmap(GL_TEXTURE_2D); // allocate the mipmaps, after texture is created.
}


/** Initialize OpenGL state. 
 *
 * Create fifteen display lists for the seven objects, and their normal vectors,
 * plus another one for the coordinate axes.
 *
 * Use Gouraud shading, and a single light source at the camera position.
 *
 * OpenGL uses a right-handed coordinate system until the perspective projection,
 * where OpenGL switches to a left-handed coordinate space for the depth buffer (z-buffer) test.
 * - The camera at the origin is looking along -Z axis in eye space, 
 * - but it is looking along +Z axis in NDC.
 * - glDepthRange is by default [0, 1] (near, far).
 * - Making the +Z axis point into the screen,
 * - and with +X to the right and +Y up. 
 * - It is a left-handed system.
 * - Changing the depth range to [1, 0] will make the system right-handed.
 *
 * @see https://www.youtube.com/watch?v=Ck1SH7oYRFM
 * @see https://webserver2.tecgraf.puc-rio.br/ftp_pub/lfm/OpenGL_Transformation.pdf
 * @see https://www.ics.uci.edu/~gopi/CS211B/opengl_programming_guide_8th_edition.pdf
 * @see http://www.songho.ca/opengl/gl_projectionmatrix.html
 *
 */
void init(void) {
   glEnable(GL_DEPTH_TEST);
   // Depth value used when the depth buffer is cleared. 
   glClearDepth(1.0);
   // Gouraud shading.
   glShadeModel (GL_SMOOTH);
   glClearColor(0.0, 0.0, 0.0, 0.0);
   //glEnable(GL_NORMALIZE);
   glEnable(GL_RESCALE_NORMAL);
   // Disable clipping against near and far planes.
   glDisable(GL_DEPTH_CLAMP);
 
   // Create light components
   GLfloat global_ambient[] = { 0.1, 0.1, 0.1, 1.0 };
   GLfloat ambient[] = { 0.05, 0.05, 0.05 };
   // diffuse light color
   GLfloat diffuseLight[] = { 0.8, 0.8, 0.8, 1.0 };
   // specular light color
   GLfloat specular[] = {1.0, 1.0, 1.0 , 1.0};

   // Assign created components to GL_LIGHT0
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
   glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
   glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
   // It's a 4D vector. If w=1, then it's a positional light source. 
   // If w=0, it's a directional light source. 
   GLfloat position[] = {0.0,0.0,zobs,1.0};
   glLightfv(GL_LIGHT0, GL_POSITION, position);

   // mcolor will be applied to both ambient and diffuse components of the material.
   // This is done for convenience because in most cases Ambient and Diffuse properties
   // of a material should be set to the same color.
   GLfloat mcolor[] = { 1.0, 0.0, 0.0, 1.0 };
   glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mcolor);

   glEnable (GL_LIGHT0);
   glEnable (GL_LIGHTING);
   // enable color tracking (specify material properties by merely calling the glColor)
   glEnable (GL_COLOR_MATERIAL);
   // set material properties which will be assigned by glColor
   glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
   // specifies the orientation of front-facing polygons. 
   // GL_CW and GL_CCW are accepted. The initial value is GL_CCW.
   glFrontFace(GL_CCW);
   // specifies the function used to compare each incoming pixel depth value 
   // with the depth value present in the depth buffer.
   glDepthFunc(GL_LEQUAL);
   glLineWidth(2.0);
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

   initTexture();

   // create 2*nObjects+1 display lists
   objects = glGenLists (theAxes+1);

   glNewList(objects, GL_COMPILE);
        torus(R1, R2, nR2, nR1, 0);
   glEndList();
   glNewList(objects+1, GL_COMPILE);
        cylinder ( 2*R1, 2*R2, nR2, nR1 , 0 );
   glEndList();
   glNewList(objects+2, GL_COMPILE);
        cone ( 2*R1, 2*R2, nR2, nR1 , 0 );
   glEndList();
   glNewList(objects+3, GL_COMPILE);
        sphere ( 1.5*R1, 1.5*R2, 1.5*R3, nR2, nR1 , 0 );
   glEndList();
   glNewList(objects+4, GL_COMPILE);
        paraboloid ( 2*R1, 2*R2, nR2, nR1 , 0 );
   glEndList();
   glNewList(objects+5, GL_COMPILE);
        hyperbolic_paraboloid ( 2*R1, 2*R1, nR1/2, nR1/2, 0 );
   glEndList();
   glNewList(objects+6, GL_COMPILE);
        hyperboloid ( R1, R2, R3, nR1/3, nR2, 0 );
   glEndList();
   glNewList(objects+7, GL_COMPILE);
        torus(R1, R2, nR2, nR1, 1);
   glEndList();
   glNewList(objects+8, GL_COMPILE);
        cylinder ( 2*R1, 2*R2, nR2, nR1, 1 );
   glEndList();
   glNewList(objects+9, GL_COMPILE);
        cone ( 2*R1, 2*R2, nR2, nR1, 1 );
   glEndList();
   glNewList(objects+10, GL_COMPILE);
        sphere ( 1.5*R1, 1.5*R2, 1.5*R3, nR2, nR1, 1 );
   glEndList();
   glNewList(objects+11, GL_COMPILE);
        paraboloid ( 2*R1, 2*R2, nR2, nR1, 1 );
   glEndList();
   glNewList(objects+12, GL_COMPILE);
        hyperbolic_paraboloid ( 2*R1, 2*R1, nR1/2, nR1/2, 1 );
   glEndList();
   glNewList(objects+13, GL_COMPILE);
        hyperboloid ( R1, R2, R3, nR1/3, nR2, 1 );
   glEndList();
   glNewList(objects+14, GL_COMPILE);
        drawAxes ( R1, R2 );
   glEndList();
}

/** Clear window, reset depth buffer and draw torus.
 *
 * Single buffer does not work, in general.
 *
 */
void display(void) {
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   if (useArcball) {
      // Apply arcball transformation. 
      glMatrixMode (GL_MODELVIEW);
      glLoadIdentity();
      gluLookAt(0, 0, zobs, 0, 0, 0, 0, 1, 0);
      arcball.applyRotationMatrix();
   } else {
      if ( spin ) 
        switch (axis) {
           case xAxis:
             glRotatef(rotAngle, 1.0,0.0,0.0);
             break;
           case yAxis:
             glRotatef(rotAngle,0.0,1.0,0.0);
             break;
           case zAxis:
             glRotatef(rotAngle,0.0,0.0,1.0);
             break;
           default:
             ;
        }
   }

   if ( showAxes ) {
      glDisable (GL_LIGHTING);
      glCallList(objects+theAxes);
      glEnable (GL_LIGHTING);
   }
   if ( showBox ) {
      glDisable (GL_LIGHTING);
      drawBox(objBoxes[OBJECT]);
      glEnable (GL_LIGHTING);
   }
   if ( drawWire ) {
      glColor3d (0.8, 0.8, 0.0);
      glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
   } else {
      glColor3d (1.0, 1.0, 1.0);
      glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      if ( applyTexture ) {
           glEnable(GL_TEXTURE_2D);
      }
   }

   glCallList(objects+OBJECT);
   glDisable(GL_TEXTURE_2D);

   glColor3d (1.0, 1.0, 1.0);
   if (drawNormals) glCallList(objects+OBJECT+nObjects);

   glutSwapBuffers();
}

/** Return the current Coordinated Universal Time (UTC) in miliseconds.
 *
 *  @return time since the Epoch (00:00:00 UTC, January 1, 1970) measured in miliseconds.
 */
int dateNow() {
   struct timeval tp;
   gettimeofday(&tp, NULL);
   int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
   return ms;
}

/**  Timer callback function to be triggered in at least delay milliseconds. 
 *
 *   The number of milliseconds is a lower bound on the time before the callback is generated. 
 *   GLUT attempts to deliver the timer callback as soon as possible after the expiration of the callback's time interval.
 *
 *   There is no support for canceling a registered callback. Instead, ignore a callback based on 
 *   its value parameter when it is triggered.
 *
 *   @param value a parameter passed to this callback function.
 *   @see https://www.opengl.org/resources/libraries/glut/spec3/node64.html
 */
void _timer(int value) {
   // Previous time.
   static int previousTimeStamp = dateNow();
   // Frames per second.
   static int numberOfFramesForFPS = 0;
   // Current time.
   int currentTime;

   // How many times we've been here in one second.
   currentTime = dateNow();
   if ( fabs(currentTime - previousTimeStamp) >= 1000 ) {
       if ( numberOfFramesForFPS > 60 ) delay *= 2;
       else if ( numberOfFramesForFPS < 30 ) delay /= 2;
       printf("\rfps = %d, delay = %d", numberOfFramesForFPS, delay);
       fflush(stdout);
       numberOfFramesForFPS = 0;
       previousTimeStamp = currentTime;
   }
   
   numberOfFramesForFPS++;

   // send redisplay event
   glutPostRedisplay();

   // call this function again in delay milliseconds
   glutTimerFunc(delay, _timer, value);
}

/** Handle window resize.
 *
 *  Sets the viewport to its new size.
 *
 *  @param w viewport width.
 *  @param h viewport height.
 */
void reshape(int w, int h) {
   double d, dx, dy, dz;
   d = getBoxSize((oType)OBJECT, &dx, &dy, &dz);

   // Set the viewport to be the entire window
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
#if 0
   // set camera angle to fit object
   viewingAngle = 2.0 * atan2(d*0.8,zobs);
   viewingAngle = toDeg(viewingAngle);
#endif
   gluPerspective(viewingAngle, (GLfloat)w/h, frontPlane, farPlane);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
#if 1
   // set camera position to fit object
   zobs = d * 0.8 / tan(toRad(viewingAngle*0.5));
#endif
   gluLookAt(0, 0, zobs, 0, 0, 0, 0, 1, 0);
   arcball.setWidthHeight(w, h);
}

/** Keyboard callback interface.
 *
 *  glutPostRedisplay() marks the current window as needing to be redisplayed.
 *  The next iteration through glutMainLoop, the window's display callback 
 *  will be called to redisplay the window's normal plane. 
 *
 *  @param key key pressed.
 *  @param x coordinate of the mouse cursor.
 *  @param y coordinate of the mouse cursor.
 *  @see https://www.opengl.org/resources/libraries/glut/spec3/node49.html
 */
void keyboard(unsigned char key, int x, int y) {
   switch (key) {
   case 'a':
   case 'A':
      showAxes = !showAxes;
      glutPostRedisplay();
      break;
   case 'b':
   case 'B':
      showBox = !showBox;
      glutPostRedisplay();
      break;
   case 'x':
   case 'X':
      if (!spin) { 
         glRotatef(rotAngle,1.0,0.0,0.0);
         glutPostRedisplay();
      }
      axis = xAxis;
      break;
   case 'y':
   case 'Y':
      if (!spin) {
         glRotatef(rotAngle,0.0,1.0,0.0);
         glutPostRedisplay();
      }
      axis = yAxis;
      break;
   case 'z':
   case 'Z':
      if (!spin) {
         glRotatef(rotAngle,0.0,0.0,1.0);
         glutPostRedisplay();
      }
      axis = zAxis;
      break;
   case 'i':
   case 'I':
      glLoadIdentity();
      gluLookAt(0, 0, zobs, 0, 0, 0, 0, 1, 0);
      axis = xAxis;
      arcball.reset();
      glutPostRedisplay();
      break;
   case 'w':
   case 'W':
      drawWire = !drawWire;
      glutPostRedisplay();
      break;
   case 'n':
   case 'N':
      drawNormals = !drawNormals;
      glutPostRedisplay();
      break;
   case 'o':
   case 'O':
      ++OBJECT;
      OBJECT %= nObjects;
      goto code;
   case '1':
   case '2':
   case '3':
   case '4':
   case '5':
   case '6':
   case '7':
      OBJECT = key - '0' - 1;
      code:
      glutSetWindowTitle(objNames[OBJECT]);
      GLint viewport[4];
      glGetIntegerv(GL_VIEWPORT, viewport);
      reshape(viewport[2], viewport[3]);
      glutPostRedisplay();
      break;
   case 'r':
   case 'R':
      useArcball = !useArcball;
      if (useArcball) {
          GLMatrix model; 
          glGetFloatv(GL_MODELVIEW_MATRIX, model);
          arcball.setViewMatrix(model);
      }
      glutPostRedisplay();
      break;
   case 's':
   case 'S':
      spin = !spin;
      rotAngle = spin ? 2 : 30; 
      glutPostRedisplay();
      break;
   case 't':
   case 'T':
      applyTexture = !applyTexture;
      glutPostRedisplay();
      break;
   case 'p':
   case 'P':
      procImage = !procImage;
      // Delete the current texture, so it can be recreated.
      glDeleteTextures(1, &texName);
      initTexture();
      glutPostRedisplay();
      break;
   case 'h':
   case 'H':
      printf ( "a, A: Draw Axes.\n"); 
      printf ( "b, B: Draw Box.\n"); 
      printf ( "x, X: Rotate about x-axis.\n"); 
      printf ( "y, Y: Rotate about y-axis.\n"); 
      printf ( "z, Z: Rotate about z-axis.\n"); 
      printf ( "w, W: Toggle wireframe view.\n");
      printf ( "i, I: Return object to initial position.\n");
      printf ( "n, N: Toggle normal display.\n");
      printf ( "o, O: Increment object display: Torus -> Cylinder -> Cone -> Sphere -> Paraboloid -> Hyperbolic Paraboloid -> Hyperboloid.\n");
      printf ( "1,2,3,4,5,6,7: Display the ith object.\n");
      printf ( "t, T: Apply texture.\n");
      printf ( "p, P: Procedural x Image Texture.\n");
      printf ( "r, R: Euler angles x Arcball.\n");
      printf ( "s, S: Toggle spin.\n");
      printf ( "h, H: Display help.\n");
      printf ( "Esc: Exit.\n");
      break;
   case 27:
      ilDeleteImages(1, &ilImage);
      printf("\n");
      exit(0);
      break;
   }
}

/** Mouse callback interface.
 *
 *  glutPostRedisplay() marks the current window as needing to be redisplayed.
 *  The next iteration through glutMainLoop, the window's display callback
 *  will be called to redisplay the window's normal plane.
 *
 *  @param state mouse button status.
 *  @param button mouse button pressed.
 *  @param _x coordinate of the mouse cursor.
 *  @param _y coordinate of the mouse cursor.
 *  @see https://www.opengl.org/resources/libraries/glut/spec3/node50.html
 */
void mouseFunc(int state, int button, int _x , int _y) {

   if (!useArcball) return;

   if ( button == GLUT_LEFT_BUTTON )
       arcball.startRotation(_x,_y);
   else
       arcball.stopRotation();


   glutPostRedisplay();
 }

/** Mouse callback interface.
 *
 *  Set the motion callback for the current window.
 *  The motion callback for a window is called when the mouse moves within
 *  the window while one or more mouse buttons are pressed. 
 *
 *  glutPostRedisplay() marks the current window as needing to be redisplayed.
 *  The next iteration through glutMainLoop, the window's display callback
 *  will be called to redisplay the window's normal plane.
 *
 *  @param _x coordinate of the mouse cursor.
 *  @param _y coordinate of the mouse cursor.
 *  @see https://www.opengl.org/resources/libraries/glut/spec3/node51.html
 */
void mouseDrag(int _x, int _y) {
   if (!useArcball) return;
   arcball.updateRotation(_x,_y);
   glutPostRedisplay();
}

/** Main program for testing.
 *
 *  Arguments: 
 *  
 *  - argc    - number of arguments.
 *  - argv[0] - Program name.
 *  - argv[1] - r1:  Distance from the center of the tube to the center of the torus.
 *  - argv[2] - r2:  Radius of the tube.
 *  - argv[3] - r3:  Radius in z direction for an ellipsoid.
 *  - argv[4] - nr1: number of samples (divisions) on circle r1.
 *  - argv[5] - nr2: number of samples (divisions) on circle r2.
 *  - argv[6] - fname: texture file name.
 */
int main(int argc, char **argv) {
   char image[256];
   glutInitWindowSize(512, 512);
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   glutCreateWindow(argv[0]);
   memset(image, '\0', sizeof(image));
   strncpy(image,"velazquez_texture_256.jpg", sizeof(image));
   if ( argc > 5 ) {
      sscanf(argv[1], "%lf", &R1);
      sscanf(argv[2], "%lf", &R2);
      sscanf(argv[3], "%lf", &R3);
      sscanf(argv[4], "%d", &nR1);
      sscanf(argv[5], "%d", &nR2);
      if ( argc > 6 ) strncpy(image, argv[6], sizeof(image));
   } else {
      printf ( "R1 = %lf\n", R1 );
      printf ( "R2 = %lf\n", R2 );
      printf ( "R3 = %lf\n", R3 );
      printf ( "nR1 = %d\n", nR1 );
      printf ( "nR2 = %d\n", nR2 );
      printf ( "Texture = %s\n", image );
   }	   

   ilInit();
   // load the file image with DevIL
   if ( (ilImage=loadImage(image)) == -1 ) {
        printf ("Can't load image file %s using DevIL \n", image);
        return -1;
    }

   init();
   arcball = Arcball();
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboard);
   glutMouseFunc(mouseFunc);
   glutMotionFunc(mouseDrag);
   glutTimerFunc(delay, _timer, 0);
   glutDisplayFunc(display);
   display();
   glutMainLoop();
   return 0;
}
