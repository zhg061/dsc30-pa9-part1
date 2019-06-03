/*
 * NAME: Zhaoyi Guo
 * PID: A15180402
 */
import java.sql.SQLOutput;
import java.util.*;

/**
 * class that implements the graph
 */
public class Graph {

    // TODO: define a data structure to store all the vertices with fast access
    static int powerFactor = 2;
    static double powerFactor1 = 1.0/2.0;
    static Double inf = Double.POSITIVE_INFINITY;
    private HashMap<Vertex, List<Edge>> myEdgeList;
    private HashMap<String, Vertex> myVertex;
    private ArrayList<Vertex> exploredVertex;
    private CostVertex currentV;

    /**
     * Constructor for Graph that set myEdgeList and myVertex
     */
    public Graph() {

        // define myEdgeList and myVertex
        myEdgeList = new HashMap<Vertex, List<Edge>>();
        myVertex = new HashMap<String, Vertex>();

    }

    /**
     * getter for vertex
     * return: ArrayList<Vertex> result
     */
    public ArrayList getVertexes() {
        // create a new array list, and add vertex to that array list
        // using for loop
        ArrayList<Vertex> result = new ArrayList<>();
        for (Vertex vertex: myEdgeList.keySet()) {
            result.add(vertex);
        }
        return result;
    }

    /**
     * get all the edges from the myEdgeList
     * @return Collection<List<Edge>>
     */
    public Collection<List<Edge>> getEdges() {
        // return all the values from the hashmap
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

        // if v is already exist, throw exception
        if (myEdgeList.containsKey(v))
            throw new IllegalArgumentException();
        // other wise, add the v to both myEdgeList and myVertex
        myEdgeList.put(v, new ArrayList<>());
        myVertex.put(v.getName(), v);
    }

    /**
     * Gets a collection of all the vertices in the graph
     *
     * @return collection of all the vertices in the graph
     */
    public Collection<Vertex> getVertices() {

        // get all the keys from the myEgeList

        return myEdgeList.keySet();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name name of the vertex object requested
     * @return vertex object associated with the name
     */
    public Vertex getVertex(String name) {

        // get the vertex by string, using myVertex

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

        // add a new edge from u to v if myVertex contains both nameU and nameV
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

        // add edge to both direction if myVertex contains both nameU and nameV
        if (myVertex.containsKey(nameU) && myVertex.containsKey(nameV)) {
            Vertex u = myVertex.get(nameU);
            Vertex v = myVertex.get(nameV);
            // create edges from both ways
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

        // calculate the distance from u to v

        return Math.pow((Math.pow((vx - ux), powerFactor) +
                Math.pow((vy - uy), powerFactor)), powerFactor1);
    }

    /**
     * Calculates the euclidean distance for all edges in the map using the
     * computeEuclideanCost method.
     */
    public void computeAllEuclideanDistances() {

        // compute all the distances, and update it to the edges from myEdgeList
        for (Vertex key : myEdgeList.keySet()) {
            for (Edge value: myEdgeList.get(key)) {
                double distance = computeEuclideanDistance(value.getSource().getX(), value.getSource().getY(),
                        value.getTarget().getX(), value.getTarget().getY());
                // set the distance
                value.setDistance(distance);
            }

        }

    }

    /**
     * Helper method to reset all the vertices before doing graph traversal algorithms
     */
    private void resetAllVertices() {

        // changes all the prev of all the vertexes to null
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
        // get the vertexes of both start and target
        Vertex start = myVertex.get(s);
        Vertex target = myVertex.get(t);
        if (s == t) {
            return;
        }

        // Create a stack for DFS
        LinkedList<Vertex> stack = new LinkedList<>();

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

        // reset all the prev to null
        resetAllVertices();
        if (s == t) {
            return;
        }
        Vertex start = myVertex.get(s);
        Vertex target = myVertex.get(t);
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
                Vertex targetV = edge.getTarget();
                //if ( adjV is not in discoveredSet )
                //            add adjV to frontierQueue
                if (targetV.visited == false) {
                    frontierQueue.addLast(targetV);
                    targetV.prev = currentV;
                    targetV.visited = true;
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

        /**
         * define cost and vertex
         * @param cost
         * @param vertex
         */
        public CostVertex(double cost, Vertex vertex) {
            this.cost = cost;
            this.vertex = vertex;
        }

        /**
         * compare the current CostVertex and another CostVertex
         * based on their cost
         * @param o
         * @return integer
         */
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
        // create a priorityQueue that contains CostVertex
        PriorityQueue<CostVertex> vertexQ = new PriorityQueue<>();
        ArrayList<CostVertex> vertexL = new ArrayList<>();
        Vertex start = myVertex.get(s);
        Vertex target = myVertex.get(t);
        //for each vertex currentV in graph,
        // currentV->distance = Infinity, currentV->predV = 0
        for (Vertex curV : myEdgeList.keySet()) {
            if (curV == start) {
                curV.dist = 0.0;
                vertexQ.add(new CostVertex(0, curV));
//                vertexL.add(new CostVertex(0, curV));
            }
            else {
                curV.dist = inf;
//                vertexL.add(new CostVertex(inf, curV));
            }
        }
        while (!vertexQ.isEmpty()) {
            CostVertex curV = vertexQ.poll();
            if (curV.vertex == target)
                return;
            // get all the adjacent edges of the current vertex
            Vertex startV = curV.vertex;
            for(Edge edge : myEdgeList.get(startV)) {
                // get the target and distance of the edge
                Vertex adjV = edge.getTarget();
                double edgeWeight = edge.getDistance();
                double newDis = edgeWeight + curV.cost;
                // create a new cost vertex with the new distance

                if (adjV.dist > newDis && adjV.visited == false) {
                    CostVertex newCV = new CostVertex(newDis, adjV);
                    adjV.dist = newDis;
                    vertexQ.add(newCV);
                    adjV.prev = curV.vertex;

                }


                // if the target vertex has not been visited
//                if (!adjV.visited) {
//                    for (int j = 0; j < vertexL.size(); j++) {
//                        if (vertexL.get(j).vertex == adjV) {
//                            //replace the original costVertex with the new costVertex
//                            vertexQ.remove(vertexL.get(j));
//                            vertexQ.add(newCV);
//                            vertexL.remove(vertexL.get(j));
//                            vertexL.add(newCV);
//                        }
//                    }
//                    adjV.prev = curV.vertex;
//                    adjV.visited = true;
//                }
//                // if the vertex has already been visited
//                else {
//                for (int j = 0; j < vertexL.size(); j++) {
//                    // if the old distance is longger than the new distance, replace the
//                    // old costVertex with the new costVertex
//                    if (vertexL.get(j).vertex == adjV && vertexL.get(j).cost > newDis) {
//                        vertexQ.remove(vertexL.get(j));
//                        vertexQ.add(newCV);
//                        vertexL.remove(vertexL.get(j));
//                        vertexL.add(newCV);
//                    }
//                }

//                }
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
    private double hValue(String cur, String goal) {

        // calculate the h value in A*

        return 0.0;
    }

    /**
     * Find the path from vertex with name s to vertex with name t, using A*
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void AStar(String s, String t) {

        // TODO

    }

    /**
     * Returns a list of edges for a path from city s to city t.
     *
     * @param s starting city name
     * @param t ending city name
     * @return list of edges from s to t
     */
    public List<Edge> getPath(String s, String t) throws OutOfMemoryError{

        // TODO
        Stack<Edge> reversedResult = new Stack<>();
        List<Edge> result = new ArrayList<Edge>();
        Vertex target = myVertex.get(t);
        while (target.getPrev() != null) {
            Edge edge = new Edge(target.getPrev(), target,
                    computeEuclideanDistance(target.getPrev().getX(),
                            target.getPrev().getY(), target.getX(), target.getY()));
            reversedResult.push(edge);
            target = target.prev;
        }
        while (!reversedResult.empty()){
            result.add(reversedResult.pop());
        }
        return result;

    }

}
