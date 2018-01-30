#!/usr/bin/env python
# coding: UTF-8
#
## \page Package5 DerivedClasses.py - An extension of the polygon class suited to the polyhedra application.
#
## @package DerivedClasses
#
#  Adds new functionalities to the polygon class.
#
#  @author Flavia Cavalcanti
#  @since 02/02/2017
#

import sys

from geometry import* 
from graph import *
import numpy as np
import matrix, geometry

## Extends the Polygon class.
#
#  Adds several new attributes, such as:
#  - isSelected.
#  - rotation matrix.
#  - texture coordinates for each vertex.
#  - a local coordinate system using the plane of the polygon,
#    and setting the polygon normal as the z axis. 
#
class SelectablePoly(Polygon):
	## Constructor.
	#
	#  @param vlist - list of 3D points for all polygons
	#  @param vIndexes - list of vertex indexes
	#  @param num - index of a face in relation to the overall polyhedron
	#
	def __init__(self, vlist, vIndexes, num):
		points = []
		for i in vIndexes:
			points.append(vlist[i])

		Polygon.__init__(self, points)

		## Did the user select the polygon.
		self.isSelected = False
		## Indexes used to build this polygon from a vertex list "vList".
		self.vIndexes = vIndexes
		## Unique polygon id.
		self.num = num

		## A vector on the polygon plane.
		self.u = (points[1]-points[0]).normalize()
		## Another vector on the polygon plane.
		self.v = self.u.crossProd(self.normal)
		k = -self.u.dotProd(points[0])
		l = -self.v.dotProd(points[0])
		m = -self.normal.dotProd(points[0])
		## Maps from the uvn (polygon plane) coordinate system to xyz (z=0 plane)
		self.V = np.matrix([[self.u.x, self.u.y, self.u.z, k], \
				   		    [self.v.x, self.v.y, self.v.z, l], \
					        [self.normal.x, self.normal.y, self.normal.z, m], 
					        [0, 0, 0, 1]])
		## Texture coordinate list: one entry per vertex
		self.texture = None
		## Rotation matrix for this polygon, to flatten it.
		self.rot = None


	def toPoints(self,lpoints):
		""""Creates a point list from a list of coordinates."""

		vlist = []
		for p in lpoints:
			vlist.append(Point(p[0], p[1], p[2]))
				
		return vlist

	def setSelection(self, selection):
		""" Set this polygon to selected/deselected. """

		self.isSelected = selection


	## Given a polygon, return the common edge with this polygon.
	#  
	#  @param pol given polygon.
	#  @return a tuple with the two vertices of the common edge, or None. 
	#
	def getSharedEdge(self, pol): 
		""" Given two polygons, return the edge that they share. 
		    Return none, if they are not adjacent by an edge. 
		"""  

		# Both polygons must have the same number of vertices to prevent conflicts
		if (pol is None) or (self.n != pol.n):
			return None    

		if False:
			edges = ()
			for v1 in self.vIndexes: 
				for v2 in pol.vIndexes:
					if (v1 == v2):
						edges += (v1,)
						if len(edges) == 2:
							break
		else:
			# intersect the two sets of indexes
			edges = tuple(set(self.vIndexes) & set(pol.vIndexes))

		if len(edges) != 2:
			return None
	
		return edges

	## Apply a given projective matrix to the vertices of this polygon.
	#
	#  @param mat projective matrix.
	#  @return a list of points corresponding to the transformed vertices.
	#
	def transformPoints(self, mat):
		""" Given a 4x4 projective matrix, return a list containing 
            the transformed points of this polygon.
		"""

		if (mat is None):
			return 

		transformed = []
		for p in self.points:
			res = matrix.dot(mat,p.np()).tolist()[0]
			transformed.append(res)
		return transformed

##
#   Main method. Used for testing.
#	
def main(argv=None):

	#vList = [Point(1.0/3,1.0/3,1.0/3), Point(2,3,-4), Point(-1,9,-7), Point(-5,-2,8)]
	#vList = [Point(-1,-1,1), Point(-1,-1,-1), Point(1,-1,-1), Point(1,-1,1)]
	size = 4
	vList =  [Point(0,0,0), Point(size,0,0), Point(size,size,0), Point(0,size,0)]
	vIndexes = [0,1,2,3]
	poly = SelectablePoly(vList, vIndexes, 0)
	print ("Poly:\n%s" % poly)
	print ("Poly.V:\n%s" %poly.V)
	for v in vList:
		x = v.np().reshape(4,1)
		print (type(x))
		print (poly.V * x)

	print ("\nPoly = %s" % poly)
	print("Poly area = %s" % poly.area())
	print("Poly normal = %s" % poly.normal)
	transf = matrix.rotate (30, 0, 0, 1)
	vtx = poly.transformPoints(transf)
	vList = poly.toPoints(vtx)
	poly2 = SelectablePoly(vList, vIndexes, 0)
	poly2.flipNormal()
	print ("\nPoly2 = %s" % poly2)
	print("Poly2 area = %s" % poly2.area())
	print("Poly2 normal = %s" % poly2.normal)
	
if __name__=="__main__":
    sys.exit(main())
