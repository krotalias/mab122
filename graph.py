#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## @package graph
#
#  A graph data structure consists of a finite (and possibly mutable) set 
#  of vertices or nodes or points, together with a set of unordered pairs of 
#  these vertices for an undirected graph or a set of ordered pairs for a directed graph. 
#  These pairs are known as edges, arcs, or lines for an undirected graph and as arrows, 
#  directed edges, directed arcs, or directed lines for a directed graph. 
#  The vertices may be part of the graph structure, or may be external entities represented 
#  by integer indices or references.
#
#  @author Flavia Cavalcanti
#  @date 12/02/2017 
#  @see http://www.python-course.eu/graphs_python.php
#
import sys

class Graph(object):
    """ A Graph Class.
        A simple Python graph class, demonstrating the essential 
        facts and functionalities of graphs.
    """

    def __init__(self, graph_dict=None):
        """ initializes a graph object 
            If no dictionary or None is given, 
            an empty dictionary will be used
        """
        if graph_dict == None:
            graph_dict = {}
        self.__graph_dict = graph_dict

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def add_edge(self, edge):
        """ assumes that edge is of type set, tuple or list; 
            between two vertices can be multiple edges! 
        """
        edge = set(edge)
        (vertex1, vertex2) = tuple(edge)
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]

    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if sys.hexversion < 0x02070000:
                    if set([neighbour, vertex]) not in edges:
                        edges.append(set([vertex, neighbour]))
                else:
                    if {neighbour, vertex} not in edges:
                        edges.append({vertex, neighbour})

        return edges

    def get(self, key):
        return self.__graph_dict.get(key)

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res


if __name__ == "__main__":

    g = { "a" : ["d"],
          "b" : ["c", "d"],
          "c" : ["a", "b"],
          "d" : ["a", "b"],
        }


    graph = Graph(g)
    print ("Complete graph:\n%s\n" % graph)

    print("Vertices of graph:")
    print(graph.vertices())

    print("Edges of graph:")
    print(graph.edges())

    print("Add vertex:")
    graph.add_vertex("z")

    print("Vertices of graph:")
    print(graph.vertices())
 
    print("Add an edge:")
    if sys.hexversion < 0x02070000:
        graph.add_edge(set(["a","z"]))
    else:
        graph.add_edge({"a","z"})
    
    print("Vertices of graph:")
    print(graph.vertices())

    print("Edges of graph:")
    print(graph.edges())

    print('Adding an edge {"x","y"} with new vertices:')
    if sys.hexversion < 0x02070000:
        graph.add_edge(set(["x","y"]))
    else:
        graph.add_edge({"x","y"})
    print("Vertices of graph:")
    print(graph.vertices())
    print("Edges of graph:")
    print(graph.edges())
