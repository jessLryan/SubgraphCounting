package subgraphembedding;

import java.util.concurrent.CopyOnWriteArrayList;

/**
 *
 * @author jessica
 */
public class Vertex {
    private CopyOnWriteArrayList<Vertex> neighbours;
    private boolean delete = false;
    
    /**
     *
     */
    public Vertex() {
        this.neighbours = new CopyOnWriteArrayList<>();
    }
    
    /**
     *
     * @param v
     */
    public void makeAdjacent(Vertex v) {
        neighbours.add(v);
    }
    
    /**
     *
     */
    public void deleteNeighbours() {
        for (Vertex v:neighbours) {
            v.delAdjacency(this);
            this.delAdjacency(v);
        }
    }
    
    /**
     *
     * @param value
     */
    public void setDeleteFlag (boolean value) {
        delete = value;
    }
    
    /**
     *
     * @return
     */
    public boolean getDeleteFlag () {
        return delete;
    }
    
    /**
     *
     * @param v
     * @return
     */
    public boolean isAdjacent(Vertex v) {
        return neighbours.indexOf(v) != -1;
    }
     
    /**
     *
     * @param v
     */
    public void delAdjacency(Vertex v) {
        neighbours.remove(v);
     }
     
    /**
     *
     * @return
     */
    public CopyOnWriteArrayList<Vertex> getNeighbours() {
        return neighbours;
    }
    
    /**
     *
     * @return
     */
    public int degree() {
        return neighbours.size();
    }
    
}