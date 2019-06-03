/*
 * NAME: Zhaoyi Guo
 * PID: A15180402
 */

import java.util.*;

/**
 * class that implements Vertex
 */
public class Vertex {

    private final String name; // the name of this vertex
    private final int x; // the x coordinates of this vertex on map
    private final int y; // the y coordinates of this vertex on map
    public double dist;
    public Vertex prev;
    public boolean visited;
    // TODO: add additional instance variables to work with different graph traversal algorithm

    private List<Edge> edgesList;

    /**
     * constructor of Vertex that sets name, x, y, and visited
     * @param name
     * @param x
     * @param y
     */
    public Vertex(String name, int x, int y) {
        // define x, y, name, and visited
        this.name = name;
        this.x = x;
        this.y = y;
        visited = false;
    }

    /**
     * get the name of the vertex
     * @return name
     */
    public String getName() {
        // get the name of the vertex
        return name;
    }

    /**
     * get the x value of the vertex
     * @return x
     */
    public int getX() {
        // get the x of the vertex
        return x;
    }

    /**
     * get the y value of the vertex
     * @return y
     */
    public int getY() {
        // get the y of the vertex
        return y;
    }

    // TODO: add necessary getters and setters for ALL your instance variable

    /**
     * get the previous vertex of the current vertex
     * @return previous
     */
    public Vertex getPrev() {
        // get the prev of the vertex
        return prev;
    }

    /**
     * set the previous vertex to vertex X
     * @param x
     */
    public void setPrev(Vertex x) {
        // set the prev to x
        prev = x;
    }

    /**
     * get the list of edges relate to the current vertex
     * @return edgeList
     */
    public List<Edge> getVertices() {
        return edgesList;
    }

    /**
     * add new adge to the current vertex
     * @param edge
     */
    public void setVertices(Edge edge) {
        // add a new edge to the list of edges
        edgesList.add(edge);
    }
    @Override
    public int hashCode() {
        // we assume that each vertex has a unique name
        return name.hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null) {
            return false;
        }
        if (!(o instanceof Vertex)) {
            return false;
        }
        Vertex oVertex = (Vertex) o;

        return name.equals(oVertex.name) && x == oVertex.x && y == oVertex.y;
    }

    public String toString() {
        return name + " (" + x + ", " + y + ")";
    }

}