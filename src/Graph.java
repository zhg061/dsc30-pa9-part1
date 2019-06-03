
import java.sql.SQLOutput;
import java.util.*;

/**
 *
 */
public class Graph {

    // TODO: define a data structure to store all the vertices with fast access
    static int powerFactor = 2;
    static double powerFactor1 = 1/2;
    static Double inf = Double.POSITIVE_INFINITY;
//    is this initailization correct
    private HashMap<Vertex, List<Edge>> myEdgeList;
    private HashMap<String, Vertex> myVertex;
    private ArrayList<Vertex> exploredVertex;

    /**
     * Constructor for Graph
     */
    public Graph() {

        // TODO
        myEdgeList = new HashMap<Vertex, List<Edge>>();
        myVertex = new HashMap<String, Vertex>();

    }

    /**
     * getter for vertex
     */
    public ArrayList getVertexes() {
        ArrayList<Vertex> result = new ArrayList<>();
        for (Vertex vertex: myEdgeList.keySet()) {
            result.add(vertex);
        }
        return result;
    }

    public Collection<List<Edge>> getEdges() {
        return myEdgeList.values();
    }

    /**
     * Adds a vertex to the graph. Throws IllegalArgumentException if given vertex
     * already exist in the graph.
     *
     * @param v vertex to be added to the graph
     * @throws IllegalArgumentException if two vertices with the same name are added.
     */
    public void addVertex(Vertex v) throws IllegalArgumentException {

        // TODO
        if (myEdgeList.containsKey(v))
            throw new IllegalArgumentException();
        myEdgeList.put(v, new ArrayList<>());
        myVertex.put(v.getName(), v);
    }

    /**
     * Gets a collection of all the vertices in the graph
     *
     * @return collection of all the vertices in the graph
     */
    public Collection<Vertex> getVertices() {

        // TODO

        return myEdgeList.keySet();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name name of the vertex object requested
     * @return vertex object associated with the name
     */
    public Vertex getVertex(String name) {

        // TODO

        return myVertex.get(name);
    }

    /**
     * Adds a directed edge from vertex u to vertex v, Throws IllegalArgumentException if one of
     * the vertex does not exist
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight weight of the edge between vertex u and v
     * @throws IllegalArgumentException if one of the vertex does not exist
     */
    public void addEdge(String nameU, String nameV, Double weight) throws IllegalArgumentException {

        // TODO
        if (myVertex.containsKey(nameU) && myVertex.containsKey(nameV)) {
            Vertex u = myVertex.get(nameU);
            Vertex v = myVertex.get(nameV);
            Edge newEdge = new Edge(u, v, weight);
            myEdgeList.get(u).add(newEdge);
        }
        else
            throw new IllegalArgumentException();

    }

    /**
     * Adds an undirected edge between vertex u and vertex v by adding a directed
     * edge from u to v, then a directed edge from v to u
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight  weight of the edge between vertex u and v
     */
    public void addUndirectedEdge(String nameU, String nameV, double weight) {

        // TODO
        if (myVertex.containsKey(nameU) && myVertex.containsKey(nameV)) {
            Vertex u = myVertex.get(nameU);
            Vertex v = myVertex.get(nameV);
            Edge newEdgeU = new Edge(u, v, weight);
            Edge newEdgeV = new Edge(v, u, weight);
            myEdgeList.get(u).add(newEdgeU);
            myEdgeList.get(v).add(newEdgeV);

        }
        else
            throw new IllegalArgumentException();

    }

    /**
     * Computes the euclidean distance between two points as described by their
     * coordinates
     *
     * @param ux (double) x coordinate of point u
     * @param uy (double) y coordinate of point u
     * @param vx (double) x coordinate of point v
     * @param vy (double) y coordinate of point v
     * @return (double) distance between the two points
     */
    public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {

        // TODO

        return Math.pow((Math.pow((vx - ux), powerFactor) +
                Math.pow((vy - uy), powerFactor)), powerFactor1);
    }

    /**
     * Calculates the euclidean distance for all edges in the map using the
     * computeEuclideanCost method.
     */
    public void computeAllEuclideanDistances() {

        // TODO
        for (Vertex key : myEdgeList.keySet()) {
            for (Edge value: myEdgeList.get(key)) {
                double distance = computeEuclideanDistance(value.getSource().getX(), value.getSource().getY(),
                        value.getTarget().getX(), value.getTarget().getY());
                value.setDistance(distance);
            }

        }

    }

    /**
     * Helper method to reset all the vertices before doing graph traversal algorithms
     */
    private void resetAllVertices() {

        // TODO
        for (Vertex curV: myEdgeList.keySet()) {
            curV.setPrev(null);
        }

    }

    /**
     * Find the path from vertex with name s to vertex with name t, using DFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void DFS(String s, String t) {
        resetAllVertices();
        // TODO
        Vertex start = myVertex.get(s);
        Vertex target = myVertex.get(t);
        if (s == t) {
            return;
        }

        // Create a stack for DFS
        LinkedList<Vertex> stack = new LinkedList<>();
//        Stack<Vertex> visitedStack = new Stack<>();

        // Push the current source node
        stack.add(start);
        while (!stack.isEmpty()) {
            Vertex currentV = stack.poll();
            currentV.visited = true;
            if (currentV == target) {
                return;
            }

            //for each vertex adjV adjacent to currentV
            //Push adjV to stack
            for (Edge value: myEdgeList.get(currentV)) {
                Vertex curVertex = value.getTarget();
                if (!curVertex.visited) {
                    curVertex.visited = true;
                    curVertex.prev = currentV;
                    stack.addLast(curVertex);
                }
            }
        }


    }


    /**
     * Find the path from vertex with name s to vertex with name t, using BFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void BFS(String s, String t) {

        // TODO
        resetAllVertices();
        if (s == t) {
            return;
        }
        Vertex start = myVertex.get(s);
        Vertex target = myVertex.get(t);
        System.out.println("current target: " + target);
        LinkedList<Vertex> frontierQueue = new LinkedList<>();
        //Push startV to frontierQueue
        frontierQueue.add(start);
        //Add  startV to discoveredSet
        while (!frontierQueue.isEmpty()) {
            //currentV = Pop frontierQueue
            Vertex currentV = frontierQueue.poll();
            currentV.visited = true;
            if (currentV == target) {
                return;
            }
            //"Visit" currentV
            for (Edge edge: myEdgeList.get(currentV)) {
                System.out.println("currentV's edges: " + edge.getTarget());
                Vertex targetV = edge.getTarget();
                //if ( adjV is not in discoveredSet )
                //            Push adjV to frontierQueue
                //            Add  adjV to discoveredSet
                if (targetV.visited == false) {
                    frontierQueue.addLast(targetV);
                    targetV.prev = currentV;
                    targetV.visited = true;

                    System.out.println(targetV + " previous vertex is: " + targetV.getPrev());
                }
            }
        }
    }

    /**
     * Helper class for Dijkstra and A*, used in priority queue
     */
    private class CostVertex implements Comparable<CostVertex> {
        double cost;
        Vertex vertex;

        public CostVertex(double cost, Vertex vertex) {
            this.cost = cost;
            this.vertex = vertex;
        }

        public int compareTo(CostVertex o) {
            return Double.compare(cost, o.cost);
        }


    }

    /**
     * Find the shortest path from vertex with name s to vertex with name t, using Dijkstra
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void Dijkstra(String s, String t) {
        resetAllVertices();
        if (s == t)
            return;
        PriorityQueue<CostVertex> unvisited = new PriorityQueue<>();
        //for each vertex currentV in graph, currentV->distance = Infinity, currentV->predV = 0
        Vertex start = myVertex.get(s);
        Vertex target = myVertex.get(t);
        for (Vertex curV: myEdgeList.keySet()) {
            if (curV == start)
                unvisited.add(new CostVertex(0, curV));
            else
                unvisited.add(new CostVertex(inf, curV));
        }
        while (!unvisited.isEmpty()) {
            CostVertex currentV = unvisited.poll();
            currentV.vertex.visited = true;
            if (currentV.vertex == target)
                return;
            for (Edge value: myEdgeList.get(currentV.vertex)) {
                Vertex adjV = value.getTarget();
                if (!adjV.visited) {
                    double edgeWeight = value.getDistance();
                    double alternativePathDistance = edgeWeight + currentV.cost;
                    CostVertex newCV = new CostVertex(alternativePathDistance, adjV);
                    unvisited.add(newCV);
                    adjV.prev = currentV.vertex;
                    adjV.visited = true;
                }
            }
        }
    }

    /**
     * Helper method to calculate the h value in A*
     *
     * @param cur the current vertex being explored
     * @param goal the goal vertex to reach
     * @return the h value of cur and goal vertices
     */
//    private double hValue(String cur, String goal) {
//
//        // TODO
//
//        return 0.0;
//    }
//
//    /**
//     * Find the path from vertex with name s to vertex with name t, using A*
//     *
//     * @param s the name of starting vertex
//     * @param t the name of targeting vertex
//     */
//    public void AStar(String s, String t) {
//
//        // TODO
//
//    }

    /**
     * Returns a list of edges for a path from city s to city t.
     *
     * @param s starting city name
     * @param t ending city name
     * @return list of edges from s to t
     */
    public List<Edge> getPath(String s, String t) throws OutOfMemoryError{

        // TODO
        List<Edge> reversedResult = new ArrayList<Edge>(myEdgeList.size());
        List<Edge> result = new ArrayList<Edge>(myEdgeList.size());
        Vertex start = myVertex.get(s);
        Vertex target = myVertex.get(t);
        while (target.getPrev() != null) {
            String x = target.getPrev().getName();
            Edge edge = new Edge(target.getPrev(), target,
                    computeEuclideanDistance(target.getPrev().getX(),
                            target.getPrev().getY(), target.getX(), target.getY()));
            reversedResult.add(edge);
            target = target.prev;
        }
        for (int i = 0; i < reversedResult.size(); i++) {
            result.add(reversedResult.get(i));
        }
        return result;

    }

}
