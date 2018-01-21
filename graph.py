#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## \page Package4 graph.py - A very simple class for creating graphs, using a single dictionary.
#
## @package graph
#
#  A graph data structure consists of a finite (and possibly mutable) set 
#  of vertices or nodes or points, together with a set of unordered pairs of 
#  these vertices for an undirected graph, or a set of ordered pairs for a directed graph. 
#  These pairs are known as edges, arcs, or lines for an undirected graph and as arrows, 
#  directed edges, directed arcs, or directed lines for a directed graph. 
#  The vertices may be part of the graph structure, or may be external entities represented 
#  by integer indices or references.
#
#  @author Flavia Cavalcanti
#  @date 12/02/2017 
#  @see http://www.python-course.eu/graphs_python.php
#  @see https://en.wikipedia.org/wiki/Graph_theory
#
import sys

class Graph(object):
    """ A Graph Class.
        A simple Python graph class, demonstrating the essential 
        facts and functionalities of graphs.
    """

    def __init__(self, graph_dict=None):
        """Initializes a graph object.
           If no dictionary or None is given, 
           an empty dictionary will be used.
        """
        if graph_dict == None:
            graph_dict = {}
        ## A dictionary holding the graph. 
        ## Vertices are the keys and the data are the set of nodes connected by an edge to a vertex.
        self.__graph_dict = graph_dict

    def vertices(self):
        """Returns the vertices of a graph."""
        return list(self.__graph_dict.keys())

    def edges(self):
        """Returns the edges of a graph."""
        return self.__generate_edges()

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def add_edge(self, edge, undirected=True):
        """Assumes that edge (a,b) is of type set, tuple or list; 
           between two vertices there can be multiple edges! 
           If the graph is undirected, two edges will be created: a-b and b-a.
        """
        edge = set(edge)
        (vertex1, vertex2) = tuple(edge)

        # if the vertex already there, it will be ignored
        self.add_vertex(vertex1)
        self.add_vertex(vertex2)
        self.__graph_dict[vertex1].append(vertex2)

        if undirected:
           self.__graph_dict[vertex2].append(vertex1)

    def getEdges(self,vertex):
        """Return the set of edges connected to a given vertex."""

        edges = []
        for neighbour in self.__graph_dict[vertex]:
            edges.append(set([vertex, neighbour]))
            # this new syntax does not work with python < 2.7
            # edges.append({vertex, neighbour})
    
        return edges   

    def __generate_edges(self):
        """A static method generating the edges of this
           graph. Edges are represented as sets 
           with one (a loop back to the vertex) or two 
           vertices.
        """
        edges = []
        for vertex in self.__graph_dict:
            for e in self.getEdges(vertex):
                if e not in edges:
                   edges.append(e)
                
        return edges

    def get(self, key):
        """Return the set of nodes (vertices) connected to the key vertex."""

        return self.__graph_dict.get(key)

    def __getitem__(self, key):
        """Indexing operator. Same as get."""
        return self.get(key)

    def __str__(self):
        """Return a string representation of this graph."""

        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res

    def __repr__(self):
        """Return the string representation of the private __graph_dict."""

        return str(self.__graph_dict)

def main():
    """"Main program for testing."""

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
    # edge can be passed as a list, set or tuple
    graph.add_edge(["a","z"])
    
    print("Vertices of graph:")
    print(graph.vertices())

    print("Edges of graph:")
    print(graph.edges())

    print("Edges connected to ""a"": %s -> %s" % (graph.getEdges("a"), graph["a"]))
    print("Edges connected to ""b"": %s -> %s" % (graph.getEdges("b"), graph.get("b")))
    print("Edges connected to ""c"": %s -> %s" % (graph.getEdges("c"), graph.get("c")))
    print("Edges connected to ""d"": %s -> %s" % (graph.getEdges("d"), graph.get("d")))
    print("Edges connected to ""z"": %s -> %s" % (graph.getEdges("z"), graph.get("z")))

    print('Adding an edge {"x","y"} with new vertices:')
    graph.add_edge(set(["x","y"]))
    print("Vertices of graph:")
    print(graph.vertices())
    print("Edges of graph:")
    print(graph.edges())
    print("Graph: %s" % graph.__repr__)

if __name__ == "__main__":
    sys.exit(main())

