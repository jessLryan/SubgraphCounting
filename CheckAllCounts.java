/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package subgraphembedding;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import java.util.Set;

/**
 *
 * @author jessica
 */
public class CheckAllCounts {
    
    public static void main(String[] args) {
        String resultsFolders = args[0];
        String noOfRuns = args[1];
        int intNoOfRuns = Integer.parseInt(noOfRuns);
        boolean allCountsMatch = true;
        
        loop: for (int i=0;i<intNoOfRuns;i++) {
            File thisRunResultsFolder = new File(resultsFolders+"/"+resultsFolders+noOfRuns);
            for (File patternFolder : thisRunResultsFolder.listFiles()) {
                HashMap<Integer, Integer> ciaranResults = new HashMap<>();
                HashMap<Integer, Integer> jessResults = new HashMap<>();
                for (File jessAndCiaranFolder : patternFolder.listFiles()) {
                    if (jessAndCiaranFolder.getAbsolutePath().contains("Ciaran")) {
                        ciaranResults= getCiaranResults(jessAndCiaranFolder);
                    }
                    else {
                        jessResults = getJessResults(jessAndCiaranFolder);
                    }
                }
                for (Map.Entry<Integer, Integer> entry : ciaranResults.entrySet()) { 
                    int lambda = entry.getKey();
                    int ciaranResult = entry.getValue();
                    if (jessResults.containsKey(lambda)&& ciaranResult!=jessResults.get(lambda)) {
                        int jessResult = jessResults.get(lambda);
                        System.out.println("results do not match:");
                        System.out.println("run number = "+i+", pattern graph = "+patternFolder.getName()+", lambda = "+lambda);
                        allCountsMatch = false;
                        break loop;
                    }
                }
            }
        }
        

        if (allCountsMatch) {
            System.out.println("All counts match");
        }
        else {
            System.out.println("Counts don't all match");
        }
    }
    
    private static HashMap<Integer,Integer> getCiaranResults(File folder) {
        HashMap<Integer,Integer> results = new HashMap<>();
        for (File resultFile : folder.listFiles()) {
            String lambdaString = resultFile.getName();
            lambdaString.replace("_","." );
            int lambda = Integer.parseInt(lambdaString);
            BufferedReader ciaranReader = SubgraphEmbedding.getReader(resultFile.getAbsolutePath());
            boolean timeout = checkForTimeoutCiaran(ciaranReader);
            if (!timeout) {
                BufferedReader ciaranReader2 = SubgraphEmbedding.getReader(resultFile.getAbsolutePath());
                int result = getResultCiaran(ciaranReader2);
                results.put(lambda, result);
            }
        }
        return results;
    }
    
    private static HashMap<Integer,Integer> getJessResults(File folder) {
        HashMap<Integer,Integer> results = new HashMap<>();
        for (File resultFile : folder.listFiles()) {
            String lambdaString = resultFile.getName();
            lambdaString.replace("_","." );
            int lambda = Integer.parseInt(lambdaString);
            BufferedReader jessReader = SubgraphEmbedding.getReader(resultFile.getAbsolutePath());
            boolean timeout = checkForTimeoutJess(jessReader);
            if (!timeout) {
                BufferedReader jessReader2 = SubgraphEmbedding.getReader(resultFile.getAbsolutePath());
                int result = getResultJess(jessReader2);
                results.put(lambda, result);
            }
        }
        return results;
    }
    
    private static boolean checkForTimeoutJess(BufferedReader jessReader) {
        try {
            //line i of the file lists the number of vertices in the graph
            String line = jessReader.readLine();
            while (line!=null&&!line.isEmpty()) {
                if (line.contains("TIMEOUT")||line.contains("FAIL")) { 
                    return true;
                }
                else {
                    line = jessReader.readLine();
                }
            }
        } catch (IOException e) {
            System.out.print("I/O Error");
            System.exit(0);
        }
        return false;
    }
    
    private static boolean checkForTimeoutCiaran(BufferedReader ciaranReader) {
        try {
            //line i of the file lists the number of vertices in the graph
            String line = ciaranReader.readLine();
            while (line!=null&&!line.isEmpty()) {
                if (line.contains("aborted")) { 
                    return true;
                }
                else {
                    line = ciaranReader.readLine();
                }
            }
        } catch (IOException e) {
            System.out.print("I/O Error");
            System.exit(0);
        }
        return false;
    }
    
    
    private static int getResultJess(BufferedReader jessReader) {
        int jessResult = -1;
        try {
            //line i of the file lists the number of vertices in the graph
            String line = jessReader.readLine();
            while (1==1) {
                if (line.contains("count")) { 
                    jessResult = Integer.parseInt(line.substring(8 , line.length()));
                    return jessResult;
                }
                else {
                    line = jessReader.readLine();
                }
            }
        } catch (IOException e) {
            System.out.print("I/O Error");
            System.exit(0);
        }
        return jessResult;
    }
    
    private static int getResultCiaran(BufferedReader ciaranReader) {
        int ciaranResult = -1;
        try {
            //line i of the file lists the number of vertices in the graph
            String line = ciaranReader.readLine();
            while (1==1) {
                if (line.contains("solution_count")) { 
                    ciaranResult = Integer.parseInt(line.substring(17 , line.length()));
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
    
}
