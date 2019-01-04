package subgraphembedding;

import java.util.List;
import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

/**
 *
 * @author jessica
 */
public class Graph {
    private ArrayList<Vertex> vertices;
        
    /**
     *
     */
    public Graph() {
        this.vertices = new ArrayList<>();
    }
    
    /**
     *
     * @param n
     */
    public Graph(int n) {
        this.vertices = new ArrayList<>();
        for (int i=0;i<n;i++) {
            Vertex v = new Vertex();
            vertices.add(v);
            for (int j=0;j<i;j++) {
                Vertex u = vertices.get(j);
                u.makeAdjacent(v);
                v.makeAdjacent(u);
            }
        }
    }
    
    public Graph(ArrayList<Vertex> inputVertices) {
        this.vertices = inputVertices;
    }
    
    /**
     *
     * @return
     */
    public int order() {
        return vertices.size();
    }
    
    /**
     *
     * @return
     */
    public int edgeCount() {
        int count = 0;
        for (Vertex v : vertices) {
            count+=v.degree();
        }
        return count/2;       
    }
    
    /**
     *
     * @param v
     * @return
     */
    public boolean contains(Vertex v) {
        return vertices.contains(v);
    }
    
    /**
     *
     * @param G2
     */
    public void addGraph(Graph G2) {
        int order = order();
        for (Vertex v_2 : G2.getVertices()) {
            Vertex newv_2 = new Vertex();
            vertices.add(newv_2);
        }
        for (Vertex v_2 : G2.getVertices()) {
            for (Vertex v_2N : v_2.getNeighbours()) {
                Vertex v_1 = vertices.get(order+G2.getIndexOfVertex(v_2)-1);
                Vertex v_1N = vertices.get(order+G2.getIndexOfVertex(v_2N)-1);
                if (!v_1.isAdjacent(v_1N)) {
                    v_1.makeAdjacent(v_1N);
                    v_1N.makeAdjacent(v_1);
                }
            }
        }
    }
        
    private List<Vertex> nonNeighbours(Vertex v) {
        List<Vertex> nonNeighb = new ArrayList<>();
        for (Vertex nv : vertices) {
            if (!v.getNeighbours().contains(nv)&&!v.equals(nv)) {
                nonNeighb.add(v);
            }
        }
        return nonNeighb;
    }
        
    /**
     *
     * @return
     */
    public int maxDeg() {
        int maxDeg = 0;
        for (Vertex v : vertices) {
            if (v.degree() > maxDeg) {
                maxDeg = v.degree();
            }
        }
        return maxDeg;
        
    }
    
    //remove all degree 0 vertices from graph

    /**
     *
     */
    public void removeDegZero() {
        for (Iterator<Vertex> itr = vertices.iterator(); itr.hasNext();) {
            Vertex v = itr.next();
            if (v.degree()==0) {
                itr.remove();
            }
        }
    }
    
    /**
     *
     * @param n
     */
    public void attachVertices(int n) {
        for (int i=0;i< n; i++) {
            Vertex v = new Vertex();
            for (Vertex u : vertices) {
                v.makeAdjacent(u);
                u.makeAdjacent(v);
            }
            vertices.add(v);
        }
    }
    
    /**
     *
     */
    public void attachVertex() {
            Vertex v = new Vertex();
            for (Vertex u : vertices) {
                v.makeAdjacent(u);
                u.makeAdjacent(v);
            }
            vertices.add(v);
    }
    
    /**
     *
     * @param k
     * @return
     */
    public int countDegAtLeast(int k) {
        int result = 0;
        for (Vertex v : vertices) {
            if(v.degree() >= k) {
                result++;
            }
        }
        return result;
    }
    
    
    private int[][] matrix() {
        int[][] matrix = new int[vertices.size()][vertices.size()];
        for (int i=0;i<matrix.length;i++) {
            for (int j=0;j<matrix[i].length;j++) {
                if (vertices.get(i).isAdjacent(vertices.get(j))) {
                    matrix[i][j] = 1;
                }
                else {
                    matrix[i][j] = 0;
                }
            }
        }
        return matrix;
    }

    
    /**
     *
     * @return
     */
    public int[] getDegreesArray() {
            int maxDeg = 0;
            for (Vertex v: vertices) {
                int vDeg = v.degree();
                if (vDeg>maxDeg) {
                    maxDeg = vDeg;
                }
            } 
            int[] degArr = new int[maxDeg+1];
            for (Vertex v: vertices) {
                int vDeg = v.degree();
                degArr[vDeg] ++;
            }
            return degArr;
        }
    
    /**
     *
     * @param filepath
     * @param folder
     * @param filenames
     */
    public void writeToFile(String filepath,String folder, String filenames) {
        String directoryName = filepath.concat(folder);
        String fileName = filenames+ ".txt";
        File directory = new File(directoryName);
        if (!directory.exists()){
            directory.mkdirs();
        }
        try{
            FileWriter fw = new FileWriter(directoryName + "/" + fileName);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write(Integer.toString(vertices.size()));
            bw.newLine();
			for (Vertex v : vertices) {
                            bw.write(Integer.toString(v.getNeighbours().size()));
                            for (Vertex u : v.getNeighbours()) {
                                bw.write(" "+Integer.toString(vertices.indexOf(u)));
                            }
                            bw.newLine();
                        }
                        bw.flush();
                        bw.close();
		} catch (IOException e){
                    System.exit(-1);
                }        
    }
    
    /**
     *
     * @param keepInd
     * @param deleteInd
     */
    public void mergeVertices(Integer keepInd,Integer deleteInd) {
        Vertex keep = vertices.get(keepInd);
        Vertex delete = vertices.get(deleteInd);
        List<Vertex> keepNeighbours = keep.getNeighbours();
        List<Vertex> deleteNeighbours = delete.getNeighbours();
        for (Vertex deleteNeighbour : deleteNeighbours) {
            if (!keepNeighbours.contains(deleteNeighbour)) {
                keep.makeAdjacent(deleteNeighbour);
                deleteNeighbour.makeAdjacent(keep);
            }
        }
    }
    
    /**
     *
     * @return
     */
    public Graph copy() {
        Graph gClone = new Graph();
        for (int i=0;i<order();i++) {
            Vertex vClone = new Vertex();
            gClone.addVertex(vClone);
        }
        for (Vertex v : vertices) {
            int vInd = vertices.indexOf(v);
            Vertex vClone = gClone.getVertexAtIndex(vInd);
            for (Vertex vNeigh : v.getNeighbours()) {
                int vNeighInd = vertices.indexOf(vNeigh);
                if (!vClone.isAdjacent(gClone.getVertexAtIndex(vNeighInd))) {
                    vClone.makeAdjacent(gClone.getVertexAtIndex(vNeighInd));
                }
                if (!gClone.getVertexAtIndex(vNeighInd).isAdjacent(vClone)) {
                    gClone.getVertexAtIndex(vNeighInd).makeAdjacent(vClone);
                }
            }
        }
        return gClone;
    }
    
    /**
     *
     * @param k
     * @return
     */
    public ArrayList<Vertex> getkHighestDeg(int k) {
        ArrayList<Vertex> kVertices = new ArrayList<>();
        int indexOfNext = 0;
        int currentDeg = 0;
        for (int i=0;i<k;i++) {
            for (Vertex v : vertices) {
                if (!kVertices.contains(v)&&v.degree()>currentDeg) {
                    currentDeg = v.degree();
                    indexOfNext = vertices.indexOf(v);
                }
            }
            kVertices.add(vertices.get(indexOfNext));
            currentDeg = 0;
        }
        return kVertices;
    }
    
    /**
     *
     * @param k
     * @return
     */
    public ArrayList<Integer> getHighestkIndices(int k) {
        ArrayList<Integer> indices = new ArrayList<>();
        int indexOfNext = 0;
        int currentDeg = 0;
        for (int i=0;i<k;i++) {
            for (Vertex v : vertices) {
                if (!indices.contains(vertices.indexOf(v))&&v.degree()>currentDeg) {
                    currentDeg = v.degree();
                    indexOfNext = vertices.indexOf(v);
                }
            }
            indices.add(indexOfNext);
            currentDeg = 0;
        }
        return indices;
    }
    
    /**
     *
     * @param k
     * @return
     */
    public int maxDegRemaining(int k) {
        ArrayList<Vertex> kVertices = getkHighestDeg(k);
        int maxDeg = 0;
        for (Vertex v : vertices) {
            if (!kVertices.contains(v)&&v.degree()>maxDeg) {
                maxDeg = v.degree();
            } 
        }
        return maxDeg;
    }
    
    /**
     *
     * @param k
     * @return
     */
    public ArrayList<Vertex> getVerticesAboveThreshold(int k) {
        ArrayList<Vertex> set = new ArrayList<>();
        for (Vertex v : vertices) {
            if(v.degree() >= k) {
                set.add(v);
            }
        }
        return set;
    }
    
    /**
     *
     * @param k
     * @return
     */
    public ArrayList<Integer> getIndicesOfVerticesAboveThreshold(int k) {
        ArrayList<Integer> set = new ArrayList<>();
        for (Vertex v : vertices) {
            if(v.degree() >= k) {
                set.add(vertices.indexOf(v));
            }
        }
        return set;
    }    
    
    /**
     *
     * @param v
     */
    public void addVertex(Vertex v) {
        vertices.add(v);
    }
    
    /**
     *
     * @param newVertices
     */
    public void addVertices(List<Vertex> newVertices) {
        for (Vertex v : newVertices) {
            vertices.add(v);
        }
    }
    
    /**
     *
     * @param oldVertices
     */
    public void addNewVertices(List<Vertex> oldVertices) {
        int oldOrder = order();
        for (Vertex v : oldVertices) {
            Vertex newv = new Vertex();
            vertices.add(newv);
        }
        for (Vertex v : oldVertices) {
            int indexV = oldVertices.indexOf(v);
            for (Vertex n:v.getNeighbours()) {
                int indexN = oldVertices.indexOf(n);
                Vertex v1 = getVertexAtIndex(oldOrder+indexV);
                Vertex v2 = getVertexAtIndex(oldOrder+indexN);
                if (!v1.isAdjacent(v2)) {
                    v1.makeAdjacent(v2);
                    v2.makeAdjacent(v1);
                }
            }
        }
    }
    
    /**
     *
     */
    public void clearVertices() {
        vertices.clear();
    }
    
    /**
     *
     * @param v
     */
    public void deleteVertex(Vertex v) {
        for (Vertex u : v.getNeighbours()) {
            u.delAdjacency(v);
        }
        vertices.remove(v);
    }
    
    /**
     *
     * @param v
     */
    public void moveVertexToFront(Vertex v) {
        vertices.remove(v);
        vertices.add(v);
    }
    
    /**
     *
     * @param v
     * @return
     */
    public int getIndexOfVertex(Vertex v) {
        return vertices.indexOf(v);
    }
    
    /**
     *
     * @return
     */
    public ArrayList<Vertex> getVertices() {
        return vertices;
    }
    
    /**
     *
     * @param i
     * @return
     */
    public Vertex getVertexAtIndex(int i) {
        return vertices.get(i);    
    }
    
    /**
     *
     * @param set
     * @return
     */
    public Graph removeFromSet(List<Integer> set) {
        Graph G = new Graph();
        ArrayList<Vertex> verticesCopy = new ArrayList<>();
        verticesCopy.addAll(vertices);
        for (int i : set) {
            Vertex v = vertices.get(i);
            verticesCopy.remove(v);
            for (Vertex u:v.getNeighbours()) {
                u.delAdjacency(v);
            }
        }
        Graph newG = new Graph(verticesCopy);
        return newG;
    }
    
    /**
     *
     * @param i
     */
    public void deleteVertexAtIndex(int i) {
        Vertex d = vertices.get(i);
        for (Vertex v : vertices) {
            if (v.getNeighbours().contains(d)) {
                v.getNeighbours().remove(d);
            }
        }
        vertices.remove(d);
    }
    
    /**
     *
     * @return
     */
    public boolean isConnected() {
        int order = vertices.size();
        if (order==0) {
            return false;
        }
        Vertex v = vertices.get(0);
        List<Vertex> reached = new ArrayList<>();
        findReachable(reached,v);
        return reached.size() == vertices.size();
    }
    
    /**
     *
     * @param M
     * @return
     */
    public List<Component> getComponents(List<List<Integer>> M) {
        List<Component> components = new ArrayList<>();
        List<Vertex> verticesCopy = new ArrayList<>();
        verticesCopy.addAll(vertices);
        while (verticesCopy.size()>0) {
            Component C = new Component();
            Vertex v = verticesCopy.get(0);
            findReachableComponents(verticesCopy,C,M,v);
            C.orderVertices();
            components.add(C);
        }
        return components;
    }
    
    //used by BFS method employed by isConnected() to add neighbours of 'current'
    //vertex to the reachable set of vertices
    private void findReachableComponents(List<Vertex> verticesCopy, Component C, List<List<Integer>> M, Vertex v) {
        Graph G = C.getGraph();
        int index = verticesCopy.indexOf(v);
        G.addVertex(v);
        C.addToM(M.get(index));
        List<Vertex> neighCopy = new ArrayList<>();
        neighCopy.addAll(v.getNeighbours());
        verticesCopy.remove(v);
        M.remove(index);       
        if (!neighCopy.isEmpty()) {
            for (Vertex u : neighCopy) {                
                if (verticesCopy.contains(u)) {
                    findReachableComponents(verticesCopy, C, M, u);
                }
            }
        }
    }
    
        //used by BFS method employed by isConnected() to add neighbours of 'current'
    //vertex to the reachable set of vertices
    private void findReachable(List<Vertex> reached, Vertex v) {
        reached.add(v);
        for (Vertex u : v.getNeighbours()) {
            if (reached.contains(u)==false) {
                findReachable(reached, u);
            }
        }
    }
    
    
    
    
        
}