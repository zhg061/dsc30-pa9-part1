/*
 * NAME: Zhaoyi Guo
 * PID: A15180402
 */

/**
 * class that implements edge
 */
public class Edge {

    private double distance; // the distance from source to target
    private final Vertex source; // the source vertex this edge starts from
    private final Vertex target; // the target vertex this edge ends at

    /**
     * constructor that sets vertex1 vertex2 and weight
     * @param vertex1
     * @param vertex2
     * @param weight
     */
    public Edge(Vertex vertex1, Vertex vertex2, double weight) {
        source = vertex1;
        target = vertex2;
        this.distance = weight;
    }

    /**
     * get the distance from vertex1 to vertex2
     * @return
     */
    public double getDistance() {
        return distance;
    }

    /**
     * set the distance from vertex1 to vertex2 with distance distance
     * @param distance
     */
    public void setDistance(double distance) {
        this.distance = distance;
    }

    public Vertex getSource() {
        return source;
    }

    /**
     * get the source of the edge, which is the target
     * @return
     */
    public Vertex getTarget() {
        return target;
    }

    /**
     * get the string of source and target
     * @return string source-target
     */
    public String toString() {
        return source + " - " + target;
    }
}