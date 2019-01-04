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
public class Stats {
    int GOrder;
    int GEdges;
    int GAvgDeg;
    int HOrder;
    int HEdges;
    int HAvgDeg;
    long endTime;
    long runTime;
    String HFile;
    String GFile;
    String startDateTime;
    String status = "F";
    String result = "TIMEOUT";
    String delta;
    String threshold;
}
