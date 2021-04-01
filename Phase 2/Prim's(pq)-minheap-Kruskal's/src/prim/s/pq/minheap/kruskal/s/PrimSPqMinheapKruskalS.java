/*
* Project Phase 2 of CPCS 324 
* section: DAR
*Group NO. : 3
* Members:
Lameer Shamsaldeen - 1806835
Nada Balsharaf - 1807769
Mariam Mahdi - 1825889
*instryctor: Dr.Bassma Alsulami

 */
package prim.s.pq.minheap.kruskal.s;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.PriorityQueue;
import javafx.util.Pair;
import java.util.LinkedList;
import java.util.Scanner;

public class PrimSPqMinheapKruskalS {

    ///*********************////
    // Class Graph to create graphs
    static class Graph {

        int vertices;
        LinkedList<Edge>[] adjList;
        ArrayList<Edge> allEdges = new ArrayList<>();

        public Graph() {
        }

        public Graph(int vertices) {
            this.vertices = vertices;
            adjList = new LinkedList[vertices + 1];
            for (int i = 0; i < vertices; i++) {
                adjList[i] = new LinkedList<>();
            }
        }

        public boolean addEdge(int source, int destination, int weight) { // addEdge ( 1 ,3 , 9)
            for (int i = 0; i < adjList[source].size(); i++) {
                if ((adjList[source].get(i).destination == destination)) {
                    return false;
                }
            }
            Edge edge = new Edge(source, destination, weight); // 1, 3, 9
            //forward edge
            adjList[source].addFirst(edge);
            edge = new Edge(destination, source, weight);
            //backward edge in undirected graph
            adjList[destination].addFirst(edge);

            allEdges.add(edge); //add to total edges
            return true;
        }

        public Graph make_Graph(int node, int edg) {
            Graph graph = new Graph(node);

            int s, d, W;
            s = 0;
            d = 1;
            W = (int) (Math.random() * 10) + 1; // from 1 to 10 
            graph.addEdge(s, d, W);
            for (int i = 2; i < node; i++) { // i number of node  
                s = i;// 3
                W = (int) (Math.random() * 10) + 1;
                d = (int) (Math.random() * (i - 1)) + 1;

                graph.addEdge(s, d, W);
            }

            for (int i = node; i < edg; i++) { // i number of edges ( 10 node they made before 9 edge so we 
                //complet on 9 edge for that we start from number of node 
                s = (int) (Math.random() * node);
                W = (int) (Math.random() * 10) + 1;
                d = (int) (Math.random() * node);
                /// check not in self loop 
                if (s == d) {
                    i--;
                    continue;
                }
                /// check if the edge is already before 
                if (!graph.addEdge(s, d, W)) {
                    i--;
                }
            }

            return graph;
        }

    }

    ///*********************////
    // Class Edge will needed when create graphs
    static class Edge {

        int source;
        int destination;
        int weight; // of edge between two vertices ( sourse and destination )

        public Edge(int source, int destination, int weight) {
            this.source = source;
            this.destination = destination;
            this.weight = weight;
        }
    }

    ///*********************////
    // Class HeapNode will needed when create Min heap 
    static class HeapNode {

        int vertex;
        int key;
    }

    ///*********************////
    // Class MinHeap will needed when implement the prim min heap  
    static class MinHeap {

        int capacity;
        int currentSize;
        HeapNode[] mH;
        int[] indexes; //will be used to decrease the key

        public MinHeap(int capacity) {
            this.capacity = capacity;
            mH = new HeapNode[capacity + 1];
            indexes = new int[capacity];
            mH[0] = new HeapNode();
            mH[0].key = Integer.MIN_VALUE;
            mH[0].vertex = -1;
            currentSize = 0;
        }

        public void display() {
            for (int i = 0; i <= currentSize; i++) {
                System.out.println(" " + mH[i].vertex + "   key   " + mH[i].key);
            }
            System.out.println("________________________");
        }

        public void insert(HeapNode x) {
            currentSize++;
            int idx = currentSize;
            mH[idx] = x;
            indexes[x.vertex] = idx;
            bubbleUp(idx);
        }

        public void bubbleUp(int pos) {
            int parentIdx = pos / 2;
            int currentIdx = pos;
            while (currentIdx > 0 && mH[parentIdx].key > mH[currentIdx].key) {
                HeapNode currentNode = mH[currentIdx];
                HeapNode parentNode = mH[parentIdx];

                //swap the positions
                indexes[currentNode.vertex] = parentIdx;
                indexes[parentNode.vertex] = currentIdx;
                swap(currentIdx, parentIdx);
                currentIdx = parentIdx;
                parentIdx = parentIdx / 2;
            }
        }

        public HeapNode extractMin() {
            HeapNode min = mH[1];
            HeapNode lastNode = mH[currentSize];
//            update the indexes[] and move the last node to the top
            indexes[lastNode.vertex] = 1;
            mH[1] = lastNode;
            mH[currentSize] = null;
            sinkDown(1);
            currentSize--;
            return min;
        }

        public void sinkDown(int k) {
            int smallest = k;
            int leftChildIdx = 2 * k;
            int rightChildIdx = 2 * k + 1;
            if (leftChildIdx < heapSize() && mH[smallest].key > mH[leftChildIdx].key) {
                smallest = leftChildIdx;
            }
            if (rightChildIdx < heapSize() && mH[smallest].key > mH[rightChildIdx].key) {
                smallest = rightChildIdx;
            }
            if (smallest != k) {

                HeapNode smallestNode = mH[smallest];
                HeapNode kNode = mH[k];

                //swap the positions
                indexes[smallestNode.vertex] = k;
                indexes[kNode.vertex] = smallest;
                swap(k, smallest);
                sinkDown(smallest);
            }
        }

        public void swap(int a, int b) {
            HeapNode temp = mH[a];
            mH[a] = mH[b];
            mH[b] = temp;
        }

        public boolean isEmpty() {
            return currentSize == 0;
        }

        public int heapSize() {
            return currentSize;
        }

    }

    ///*********************////
    // Class ResultSet will needed when svae MST for prim algorithms
    static class ResultSet {

        int parent;
        int weight;
    }

    ///*********************////
    // Class for application of Prim by Priority Queue  
    static class PPQ {

        static long startTime;

        // function to create minimum spanning tree
        public static void primMST(Graph g) {
            startTime = System.currentTimeMillis();
            boolean[] mst = new boolean[g.vertices];

            // to save final resulte of MST
            ResultSet[] resultSet = new ResultSet[g.vertices];
            int[] key = new int[g.vertices];  //keys used to store the key to know whether priority queue update is required

            //Initialize all the keys to infinity and
            //initialize resultSet for all the vertices
            for (int i = 0; i < g.vertices; i++) {
                key[i] = Integer.MAX_VALUE;
                resultSet[i] = new ResultSet();
            }

            //Initialize priority queue
            //override the comparator to do the sorting based keys
            PriorityQueue<Pair<Integer, Integer>> pq = new PriorityQueue<>(g.vertices, new Comparator<Pair<Integer, Integer>>() {
                @Override
                public int compare(Pair<Integer, Integer> p1, Pair<Integer, Integer> p2) {
                    //sort using key values
                    int key1 = p1.getKey();
                    int key2 = p2.getKey();
                    return key1 - key2;
                }
            });

            //create the pair for for the first index, 0 key 0 index
            key[0] = 0;
            Pair<Integer, Integer> p0 = new Pair<>(key[0], 0);
            //add it to pq
            pq.offer(p0);

            resultSet[0] = new ResultSet();
            resultSet[0].parent = -1;

            //while priority queue is not empty
            while (!pq.isEmpty()) {
                //extract the min
                Pair<Integer, Integer> extractedPair = pq.poll();

                //extracted vertex
                int extractedVertex = extractedPair.getValue();
                mst[extractedVertex] = true;

                //iterate through all the adjacent vertices and update the keys
                LinkedList<Edge> list = g.adjList[extractedVertex];
                for (int i = 0; i < list.size(); i++) {
                    Edge edge = list.get(i);
                    //only if edge destination is not present in mst
                    if (mst[edge.destination] == false) {
                        int destination = edge.destination;
                        int newKey = edge.weight;
                        //check if updated key < existing key, if yes, update if
                        if (key[destination] > newKey) {
                            //add it to the priority queue
                            Pair<Integer, Integer> p = new Pair<>(newKey, destination);
                            pq.offer(p);
                            //update the resultSet for destination vertex
                            resultSet[destination].parent = extractedVertex;
                            resultSet[destination].weight = newKey;
                            //update the key[]
                            key[destination] = newKey;
                        }
                    }
                }
            }
            //print mst 
            printMST(resultSet, g);
        }

        // called this function by primMST to print the edges of MST
        public static void printMST(ResultSet[] resultSet, Graph g) {
            long endTime = System.currentTimeMillis(); // before starting the print we must stopped the time
            int total_min_weight = 0;
            System.out.println("Minimum Spanning Tree: ");
            for (int i = 0; i < g.vertices; i++) { ///////////////////////////////////////////// can change to another ?
                System.out.println("Edge: ( " + i + " => " + resultSet[i].parent
                        + ", weight: " + resultSet[i].weight + " )");
                total_min_weight += resultSet[i].weight;
            }
            System.out.println("Total Cost: " + total_min_weight);

            // to calculate and save the time in suitable format
            NumberFormat formatter = new DecimalFormat("#0.00000");
            System.out.println("Total time in Prim's algorithm using priority-queue is " + formatter.format((endTime - startTime) / 1000d) + " seconds");
        }
    }

    ///*********************////
    // Class for application of Prim by Min Heap  
    static class PMH { //PMH => Prim Min Heap

        static long startTime;

        public static void primMST(Graph g) {
            startTime = System.currentTimeMillis();
            boolean[] inHeap = new boolean[g.vertices];
            ResultSet[] resultSet = new ResultSet[g.vertices];
            //keys[] used to store the key to know whether min hea update is required
            int[] key = new int[g.vertices];
            //initialize and create heapNode for all the vertices 
            HeapNode[] heapNodes = new HeapNode[g.vertices];
            for (int i = 0; i < g.vertices; i++) {
                heapNodes[i] = new HeapNode();
                heapNodes[i].vertex = i;
                heapNodes[i].key = Integer.MAX_VALUE;
                resultSet[i] = new ResultSet();
                resultSet[i].parent = -1;
                inHeap[i] = true;
                key[i] = Integer.MAX_VALUE;
            }

            //decrease the key for the first index
            heapNodes[0].key = 0;

            //add all the vertices to the MinHeap
            MinHeap minHeap = new MinHeap(g.vertices);
            //add all the vertices to priority queue
            for (int i = 0; i < g.vertices; i++) {
                minHeap.insert(heapNodes[i]);
            }

            //while minHeap is not empty
            while (!minHeap.isEmpty()) {
                //extract the min
                HeapNode extractedNode = minHeap.extractMin();
                //**used to know which vertex is added to vt
                System.out.println("vertex added to vt: (vertex: " + extractedNode.vertex + ", weight: " + extractedNode.key + ")");
                //extracted vertex
                int extractedVertex = extractedNode.vertex;
                inHeap[extractedVertex] = false;

                //iterate through all the adjacent vertices
                LinkedList<Edge> list = g.adjList[extractedVertex];
                for (int i = 0; i < list.size(); i++) {
                    Edge edge = list.get(i);
                    //only if edge destination is present in heap
                    if (inHeap[edge.destination]) {
                        int destination = edge.destination;
                        int newKey = edge.weight;
                        //**print the fring of the added vertex 
                        System.out.println("v-vt: ( " + edge.source + " => " + destination + ", weight: " + newKey + " )");
                        //check if updated key < existing key, if yes, update if
                        if (key[destination] > newKey) {
                            decreaseKey(minHeap, newKey, destination);
                            //update the parent node for destination
                            resultSet[destination].parent = extractedVertex;
                            resultSet[destination].weight = newKey;
                            key[destination] = newKey;
                        }
                    }
                }
            }
            //print mst
            printMST(resultSet, g);
        }

        public static void decreaseKey(MinHeap minHeap, int newKey, int vertex) {

            //get the index which key's needs a decrease;
            int index = minHeap.indexes[vertex];

            //get the node and update its value
            HeapNode node = minHeap.mH[index];
            node.key = newKey;
            minHeap.bubbleUp(index);
        }

        public static void printMST(ResultSet[] resultSet, Graph g) {
            long endTime = System.currentTimeMillis();
            int total_min_weight = 0;
            System.out.println("Minimum Spanning Tree: ");
            for (int i = 1; i < g.vertices; i++) {
                System.out.println("Edge: (" + i + " => " + resultSet[i].parent
                        + ", weight: " + resultSet[i].weight + " )");
                total_min_weight += resultSet[i].weight;
            }
            System.out.println("Total Cost: " + total_min_weight);

            NumberFormat formatter = new DecimalFormat("#0.00000");
            System.out.println("Total time in Prim's algorithm using min-heap is " + formatter.format((endTime - startTime) / 1000d) + " seconds");
        }
    }

    ///*********************////
    // Class for application of kruksal 
    static class kruksal {

        static long startTime; // to calculate the time will take.

        // function to create minimum spanning tree
        public static void kruskalMST(Graph g) {
            startTime = System.currentTimeMillis();

            // initialize  Priority Queue 
            PriorityQueue<Edge> pq = new PriorityQueue<>(g.allEdges.size(), Comparator.comparingInt(o -> o.weight));

            //add all the edges to priority queue, //sort the edges on weights
            for (int i = 0; i < g.allEdges.size(); i++) {
                pq.add(g.allEdges.get(i));
            }

            //create a parent []
            int[] parent = new int[g.vertices];

            //makeset for each vertices with itself
            makeSet(parent, g);

            // to save final resulte
            ArrayList<Edge> mst = new ArrayList<>();

            //run loop until ( vertices - 1 )edges is done
            int index = 0;
            while (index < g.vertices - 1) {
                Edge edge = pq.remove(); // return and delete from Queue 
                //check if adding this edge creates a cycle or not
                int x_set = find(parent, edge.source);
                int y_set = find(parent, edge.destination);

                if (x_set == y_set) {
                    //ignore, will create cycle
                } else {
                    //add it to our final result
                    mst.add(edge);
                    index++; // increas the loop 
                    union(parent, x_set, y_set); // union between the sets
                }
            }
            //print MST
            printGraph(mst);
        }

        // called this function by kruskalMST to print the edges of MST
        public static void printGraph(ArrayList<Edge> edgeList) {
            long endTime = System.currentTimeMillis(); // before starting the print we must stopped the time
            System.out.println("Minimum Spanning Tree:- ");
            int total_min_weight = 0;
            for (int i = 0; i < edgeList.size(); i++) {
                Edge edge = edgeList.get(i);
                System.out.println("Edge: " + i + " (source: " + edge.source
                        + ", destination: " + edge.destination
                        + ", weight: " + edge.weight + ")");
                total_min_weight += edge.weight;
            }
            System.out.println("Total minimum key: " + total_min_weight);
            // to calculate and save the time in suitable format
            NumberFormat formatter = new DecimalFormat("#0.00000");
            System.out.println("Total time in Kruskal's algorithm is " + formatter.format((endTime - startTime) / 1000d) + " seconds");
        }

        public static void makeSet(int[] parent, Graph g) {
            //Make set- creating a new element with a parent pointer to itself.
            for (int i = 0; i < g.vertices; i++) {
                parent[i] = i;
            }
        }

        public static int find(int[] parent, int vertex) {
            //chain of parent pointers from x upwards through the tree
            // until an element is reached whose parent is itself
            if (parent[vertex] != vertex) {
                return find(parent, parent[vertex]);
            };
            return vertex;
        }

        public static void union(int[] parent, int x, int y) {
            int x_set_parent = find(parent, x);
            int y_set_parent = find(parent, y);
            //make x as parent of y
            parent[y_set_parent] = x_set_parent;
        }
    }

    ///*********************////
    ///***** The Main ******//// 
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        System.out.println("----------------------------------- Phase 2 -----------------------------------");
        inputCases();
        Graph gg = new Graph();
        int choice = sc.nextInt(); // read from user the choise
        System.out.println("-------------------------------------------------------------------------------");
        resultOfCases(choice, gg); // send the choise and apply the algorithms
    }

    ///*********************////
    ///***** choise casse ******//// 
    public static void resultOfCases(int x, Graph g) {
        // choose which case need to present 

        switch (x) {
            case 1: // when n=1000 , m= 10000
                g = g.make_Graph(1000, 10000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 2: // when n=1000 , m= 15000
                g = g.make_Graph(1000, 15000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 3: // when n=1000 , m= 25000
                g = g.make_Graph(1000, 25000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 4: // when n=5000 , m= 15000
                g = g.make_Graph(5000, 15000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 5: // when n=5000 , m= 25000
                g = g.make_Graph(5000, 25000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 6: // when n=10000 , m= 15000
                g = g.make_Graph(10000, 15000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 7: // when n=10000 , m= 25000
                g = g.make_Graph(10000, 25000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 8: // when n=20000 , m= 200000
                g = g.make_Graph(20000, 200000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 9: // when n=20000 , m= 300000
                g = g.make_Graph(20000, 300000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            case 10: // when n=50000 , m= 1000000
                g = g.make_Graph(50000, 1000000);
                System.out.println("Prim's Algorithm using Min-Heap:");
                PMH.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nPrim's Algorithm using Priority Queue:");
                PPQ.primMST(g);
                System.out.println("***********************************************************************"
                        + "\nKruksal Algorithm :");
                kruksal.kruskalMST(g);
                break;
            default:
                System.out.println("Wrong choice, try again");
        }

    }

    ///*********************////
    ///***** Main Menue ******//// 
    public static void inputCases() {
        System.out.println("Choose one of the following cases (n: #vertices, m: #edges)"
                + "\nA case is implemented by Prim's Alggorithm using Min-Heap and Priority Queue"
                + "\nand implemented by Kruksal Algorithm "
                + "\n1. n=1000, m=10000"
                + "\n2. n=1000, m=15000"
                + "\n3. n=1000, m=25000"
                + "\n4. n=5000, m=15000"
                + "\n5. n=5000, m=25000"
                + "\n6. n=10000, m=15000"
                + "\n7. n=10000, m=25000"
                + "\n8. n=20000, m=200000"
                + "\n9. n=20000, m=300000"
                + "\n10. n=50000, m=1000000"
                + "\n>> Enter your choice: ");
    }

}
