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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jessica
 */
public class ExtractRuntime {
    
    public static void main(String[] args) {
        String lambda = args[0];
        String patternGraph = args[1];
        String outputFile = args[2];
        String n = args[3];
        String who = args[4];


        BufferedReader reader = SubgraphEmbedding.getReader(outputFile);
        String runtime = "";
        
        if (who.contains("Ciaran")) {
            runtime = getRuntimeCiaran(reader);
        } 
        else {
            runtime = getRuntimeJess(reader);
        }
        writeToOutputFile(lambda,patternGraph,n,who,runtime);
    }
    
    
    public static String getRuntimeJess(BufferedReader jessReader) {
        String jessRuntime = "";
        try {
            //line i of the file lists the number of vertices in the graph
            String line = jessReader.readLine();
            while (1==1) {
                if (line.contains("runtime")) { 
                    jessRuntime = line.substring(10 , line.length());
                    return jessRuntime;
                }
                else {
                    line = jessReader.readLine();
                }
            }
        } catch (IOException e) {
            System.out.print("I/O Error");
            System.exit(0);
        }
        return jessRuntime;
    }
    
    
    
    
    public static String getRuntimeCiaran(BufferedReader ciaranReader) {
        String ciaranResult = "";
        try {
            //line i of the file lists the number of vertices in the graph
            String line = ciaranReader.readLine();
            while (1==1) {
                if (line.contains("runtime")) { 
                    ciaranResult = line.substring(10 , line.length());
                    return ciaranResult;
                }
                else {
                    line = ciaranReader.readLine();
                }
            }
        } catch (IOException e) {
            System.out.print("I/O Error");
            System.exit(0);
        }
        return ciaranResult;
    }
    
    public static void writeToOutputFile(String lambda, String patternGraph, String n, String who, String runtime) {
        File file = new File("/home/jessica/Documents/PhD/SubgraphEmbedding/Runtimes/"+n+"/"+patternGraph+"/"+who);
        try {
            file.createNewFile();
        } catch (IOException ex) {
            Logger.getLogger(ExtractRuntime.class.getName()).log(Level.SEVERE, null, ex);
        }
        try{
            FileWriter fw = new FileWriter(file,true);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write(lambda+" , "+runtime);
            bw.newLine();
            bw.flush();
            bw.close();

	} catch (IOException e){
                    System.exit(-1);
        }   
    }
    
    
}
