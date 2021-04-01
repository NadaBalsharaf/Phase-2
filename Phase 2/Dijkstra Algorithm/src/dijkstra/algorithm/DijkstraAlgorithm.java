/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dijkstra.algorithm;

import java.io.*;
import javafx.util.Pair;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Scanner;

/**
 *
 *
 */
public class DijkstraAlgorithm {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException {
        // number of cities
        int vertices = 12;
        // creating graph with number of cities
        Graph graph = new Graph(vertices);
        // reading the distance between cities from file and adding them to the graph
        File file = new File("Cities.txt");
        try (Scanner scan = new Scanner(file)) {
            while (scan.hasNextLine()) {
                graph.addEdge(scan.nextInt(), scan.nextInt(), scan.nextInt());
            }
        }
        // solving start from 0 : jeddah
        graph.dijkstra_GetMinDistances(0);
    }

    static class Edge {

        int source;
        int destination;
        int weight;

        public Edge(int source, int destination, int weight) {
            this.source = source;
            this.destination = destination;
            this.weight = weight;
        }
    }

    static class Graph {

        int vertices;
        LinkedList<Edge>[] adjacencylist;

        Graph(int vertices) {
            this.vertices = vertices;
            adjacencylist = new LinkedList[vertices];
            //initialize adjacency lists for all the vertices
            for (int i = 0; i < vertices; i++) {
                adjacencylist[i] = new LinkedList<>();
            }
        }

        public void addEdge(int source, int destination, int weight) {
            Edge edge = new Edge(source, destination, weight);
            adjacencylist[source].addFirst(edge);
            //for undirected graph
            edge = new Edge(destination, source, weight);
            adjacencylist[destination].addFirst(edge);
        }

        public void dijkstra_GetMinDistances(int sourceVertex) {

            boolean[] SPT = new boolean[vertices];
            //distance used to store the distance of vertex from a source
            int[] distance = new int[vertices];

            //Initialize all the distance to infinity
            for (int i = 0; i < vertices; i++) {
                distance[i] = Integer.MAX_VALUE;
            }
            //Initialize priority queue
            //override the comparator to do the sorting based keys
            PriorityQueue<Pair<Integer, Integer>> pq = new PriorityQueue<>(vertices, new Comparator<Pair<Integer, Integer>>() {
                @Override
                public int compare(Pair<Integer, Integer> p1, Pair<Integer, Integer> p2) {
                    //sort using distance values
                    int key1 = p1.getKey();
                    int key2 = p2.getKey();
                    return key1 - key2;
                }
            });
            //create the pair for for the first index, 0 distance 0 index
            distance[0] = 0;
            Pair<Integer, Integer> p0 = new Pair<>(distance[0], 0);
            //add it to pq
            pq.offer(p0);

            //while priority queue is not empty
            while (!pq.isEmpty()) {
                //extract the min
                Pair<Integer, Integer> extractedPair = pq.poll();

                //extracted vertex
                int extractedVertex = extractedPair.getValue();
                if (SPT[extractedVertex] == false) {
                    SPT[extractedVertex] = true;

                    //iterate through all the adjacent vertices and update the keys
                    LinkedList<Edge> list = adjacencylist[extractedVertex];
                    for (int i = 0; i < list.size(); i++) {
                        Edge edge = list.get(i);
                        int destination = edge.destination;
                        //only if edge destination is not present in mst
                        if (SPT[destination] == false) {
                            ///check if distance needs an update or not
                            //means check total weight from source to vertex_V is less than
                            //the current distance value, if yes then update the distance
                            int newKey = distance[extractedVertex] + edge.weight;
                            int currentKey = distance[destination];
                            if (currentKey > newKey) {
                                Pair<Integer, Integer> p = new Pair<>(newKey, destination);
                                pq.offer(p);
                                distance[destination] = newKey;
                            }
                        }
                    }
                }
            }
            //print Shortest Path Tree
            printDijkstra(distance, sourceVertex);
        }

        public void printDijkstra(int[] distance, int sourceVertex) {
            String cities[] = {"Jeddah", "Makkah", "Madinah", "Riyadh", "Dammam", "Taif", "Abha", "Tabuk", "Qasim", "Hail", "Jizan", "Najran"};
            System.out.println("------------------ Dijkstra Algorithm -------------------- "
                    + "\n distance from vertex (0: jeddah => all cities):");
            for (int i = 0; i < vertices; i++) {
                System.out.println("(From Jeddah " + sourceVertex + " => to " + i + " " + cities[i] + ", weight: " + distance[i] + ")");
            }
        }
    }

}
