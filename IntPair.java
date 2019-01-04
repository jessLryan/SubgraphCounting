/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package subgraphembedding;

/**
 *
 * @author jessica
 */
public class IntPair {
    final int x;
    final int y;
    IntPair(int x, int y) {this.x=x;this.y=y;}
    
    public int getFirst() {
        return x;
    }
    
    public int getSecond() {
        return y;
    }
}
