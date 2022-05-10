#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## \page Package6 glpick.py - Utilities needed for picking a face.
#
## @package glpick
#
#   Face selection using raycast.
#
#   @author Flavia Cavalcanti
#   @date 01/02/2017
#
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

from geometry import *

##
#  maps a point from screen coordinates to WC.
#
#  @param pt given point.
#  @return pt in world coordinates.
#
def unProject(pt):
    modelMatrix = glGetDoublev(GL_MODELVIEW_MATRIX)
    projMatrix = glGetDoublev(GL_PROJECTION_MATRIX)
    viewport = glGetIntegerv(GL_VIEWPORT)

    wcx, wcy, wcz = gluUnProject(
        pt.x, pt.y, pt.z, modelMatrix, projMatrix, viewport)

    return Point(wcx, wcy, wcz)

## Returns a ray starting at the mouse position.
#
#  @param x coordinate of the mouse.
#  @param y coordinate of the mouse.
#  @return a Line object.
#
def pickFace(x, y):
    """Returns a ray starting at the mouse position."""

    p1 = Point(x, y, 0.0)
    p2 = Point(x, y, 1.0)

    p1 = unProject(p1)
    p2 = unProject(p2)

    return Line(p1, p2)

## Main program for testing.
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # pass arguments to init
    glutInit(argv)

    # Select type of Display mode:
    #  Double buffer
    #  RGBA color
    # Alpha components supported
    # Depth buffer
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)

    Height = 480
    Width = 640

    # get a 640 x 480 window
    glutInitWindowSize(Width, Height)

    # the window starts at the upper left corner of the screen
    glutInitWindowPosition(0, 0)

    # Okay, like the C version we retain the window id to use when closing, but for those of you new
    # to Python, remember this assignment would make the variable local and not global
    # if it weren't for the global declaration at the start of main.
    window = glutCreateWindow(b"Lesson 48: NeHe ArcBall Rotation Tutorial")

    glViewport(0, 0, Width, Height)
    gluPerspective(45.0, float(Width) / float(Height), 1, 100.0)

    l4 = pickFace(100, 100)
    print(l4)


if __name__ == '__main__':
    sys.exit(main())
