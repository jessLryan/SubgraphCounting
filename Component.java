
package subgraphembedding;

import java.util.List;
import java.util.ArrayList;

//contains a graph C and associated mappable set M corresponding to some host 
//graph G

/**
 *
 * @author jessica
 */
public class Component {
    
    private final Graph C;
    private List<List<Integer>> M;

    /**
     *
     */
    public Component() {
        this.C = new Graph();
        this.M = new ArrayList<>();
    }
    
    /**
     *
     * @param C
     * @param M
     */
    public Component(Graph C, List<List<Integer>> M) {
        this.C = C;
        this.M = M;
    }
    
    //remove vertices in C (and their corresponding mappable sets in M) at indices in 'set'

    /**
     *
     * @param set
     * @return
     */
    public Component removeFromSet(List<Integer> set) {
        Graph H = C.removeFromSet(set);
        List<List<Integer>> MCopy = SubgraphEmbedding.cloneMapSet(M);
        for (int i=0;i<set.size();i++) {
            int lessThan = 0;
            for (int j=0;j<i;j++) {
                if (set.get(j)<set.get(i)) {
                    lessThan ++;
                }
            }
            MCopy.remove(set.get(i)-lessThan);
        }
        Component CNew = new Component(H,MCopy);
        return CNew;
    }
    
    /**
     *
     * @return
     */
    public Component copy() {
        Graph copyG = C.copy();
        List<List<Integer>> copyM = new ArrayList<>();
        for (List<Integer> toCopy : M) {
            List<Integer> copyTo = new ArrayList<>();
            copyTo.addAll(toCopy);
            copyM.add(copyTo);
        }
        Component compCopy = new Component(copyG,copyM);
        return compCopy;
    }
    
    /**
     *
     * @return
     */
    public Graph getGraph() {
        return C;
    }
    
    /**
     *
     * @return
     */
    public List<List<Integer>> getM() {
        return M;
    } 
     
    /**
     *
     * @param v
     */
    public void addVertex(Vertex v) {
        C.addVertex(v);
    }
        
    /**
     *
     * @param mapSet
     */
    public void addToM(List<Integer> mapSet) {
        M.add(mapSet);
    }
    
    /**
     *
     * @param v
     * @return
     */
    public boolean contains(Vertex v) {
        return C.contains(v)==true;
    }
    
    //order vertices in C so that each vertex is preceeded by at least 
    //one of its neighbours

    /**
     *
     */
    public void orderVertices() {
        if (C.order()<2) {
            return;   
        }
        Vertex current = C.getVertexAtIndex(0);
        List<Vertex> seen = new ArrayList<>();
        List<List<Integer>> seenMap = new ArrayList<>();
        seen.add(current);
        seenMap.add(M.get(0));
        addToOrder(seen, seenMap, current);      
        C.clearVertices();
        C.addVertices(seen);
        M = seenMap;
    }
    
    //used by orderVertices() to add neighbours of 'current' vertex to 'seen'. Iteratively
    //adds all vertices in C to 'seen' in such a way that each vertex is preceeded by at least
    //one of its neighbours
    private void addToOrder(List<Vertex> seen, List<List<Integer>> seenMap, Vertex current) {
        for (Vertex next : current.getNeighbours()) {
                    if (seen.contains(next) == false) {
                        current = next;
                        seen.add(current);
                        int currentIndex = C.getIndexOfVertex(current);
                        seenMap.add(M.get(currentIndex));
                        addToOrder(seen, seenMap, current);
                        break;
                    }
            }
    }
    
}