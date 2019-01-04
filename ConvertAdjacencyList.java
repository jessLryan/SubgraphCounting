/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package subgraphembedding;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jessica
 */
public class ConvertAdjacencyList {
    public static void main(String[] args) {
        BufferedReader inputFile = SubgraphEmbedding.getReader(args[0]);
        BufferedReader inputFileToFindOrder = SubgraphEmbedding.getReader(args[0]);
        int order = getOrder(inputFileToFindOrder);
        //String filename = args[1];
        //String n = args[2];
        //readAndWrite(order, inputFile,"/home/jessica/Documents/PhD/SubgraphEmbedding/testGraphs/"+n"/Converted/",filename);
        
        readAndWrite(order, inputFile);
        
    }

    
private static int getOrder(BufferedReader reader) {
    int n = 0;
        try {
            String line = reader.readLine();
            line = reader.readLine();
            while (line != null && !line.isEmpty()) {
                String firstVertex = "";
                String secondVertex = "";
                int i=1;
                while (Character.isDigit(line.charAt(i))) {
                    firstVertex = firstVertex + line.charAt(i);
                    i++;
                }
                i++;
                while (i<line.length()&&Character.isDigit(line.charAt(i))) {
                    secondVertex = secondVertex + line.charAt(i);
                    i++;
                }
                int v1 = Integer.parseInt(firstVertex);
                int v2 = Integer.parseInt(secondVertex);
                if (v1>n) {
                    n = v1;    
                }
                if (v2>n) {
                    n = v2;
                }
                line = reader.readLine();
            }
            
        } catch (IOException ex) {
            Logger.getLogger(ConvertAdjacencyList.class.getName()).log(Level.SEVERE, null, ex);
        }
        return n+1;
}  
    
    

  private static void readAndWrite(int order, BufferedReader reader) {
//private static void readAndWrite(int order, BufferedReader reader, String directoryName, String filename) {
//        String fileName = filename;
//        File directory = new File(directoryName);
//        if (!directory.exists()){
//            directory.mkdirs();
//        }
        try{
//            FileWriter fw = new FileWriter(directoryName + "/" + fileName);
//            FileWriter fwCiaran = new FileWriter(directoryName + "/" + fileName+"Ciaran");     

            FileWriter fw = new FileWriter("/home/jessica/Documents/SubgraphEmbedding/myhost");

            BufferedWriter bw = new BufferedWriter(fw);
//            BufferedWriter bwCiaran = new BufferedWriter(fwCiaran);
            String line = reader.readLine();
            bw.write(Integer.toString(order));
            bw.newLine();
//            bwCiaran.write(Integer.toString(order));
//            bwCiaran.newLine();
            List<List<Integer>> adjacencies = new ArrayList<>();
            for (int i=0;i<order;i++) {
                List<Integer> list = new ArrayList<>();
                adjacencies.add(list);
            }
            line = reader.readLine();
            while (line != null && !line.isEmpty()) {
                String firstVertex = "";
                String secondVertex = "";
                int i=1;
                while (Character.isDigit(line.charAt(i))) {
                    firstVertex = firstVertex + line.charAt(i);
                    i++;
                }
                i++;
                while (i<line.length()&&Character.isDigit(line.charAt(i))) {
                    secondVertex = secondVertex + line.charAt(i);
                    i++;
                }
                int v1 = Integer.parseInt(firstVertex);
                int v2 = Integer.parseInt(secondVertex);
                if (v1!=v2) {
                    List<Integer> v1List = adjacencies.get(v1);
                    List<Integer> v2List = adjacencies.get(v2);
                    if (!v1List.contains(v2)) {
                        v1List.add(v2);
                    }
                    if (!v2List.contains(v1)) {
                        v2List.add(v1);
                    }
                
                }
                line = reader.readLine();
            }
            for (int j=0;j<adjacencies.size();j++) {
                List<Integer> adjList = adjacencies.get(j);
                List<Integer> adjListNew = new ArrayList<>();
                adjListNew.add(adjList.size());
                adjListNew.addAll(adjList);
                adjacencies.set(j, adjListNew);
            }
            for (int k=0;k<adjacencies.size();k++) {
                List<Integer> adjList = adjacencies.get(k);
                bw.write("newVertex");
                for (Integer neigh : adjList) {
                    bw.write(Integer.toString(neigh));
                    bw.write(" ");
//                    bwCiaran.write(Integer.toString(neigh));
//                    bwCiaran.write(" ");
                    
                }
                bw.newLine();
//                bwCiaran.newLine();
            }
            bw.flush();
            bw.close();
//            bwCiaran.flush();
//            bwCiaran.close();
	} catch (IOException e){
                    System.exit(-1);
         }   
    }     



}