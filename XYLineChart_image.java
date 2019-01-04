/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package subgraphembedding;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.util.Pair;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.ChartUtilities; 

public class XYLineChart_image {

   public static void main( String[ ] args )throws Exception {
      BufferedReader jessFile = getReader(args[0]);
      BufferedReader ciaranFile = getReader(args[1]);
      String patternGraph = args[2];
      final List<DoubleIntPair> jessTimes = extractLambdaTimePair(jessFile);
      final List<DoubleIntPair> ciaranTimes = extractLambdaTimePair(ciaranFile);
      
      final XYSeries jessSeries = new XYSeries("PC");
      final XYSeries ciaranSeries = new XYSeries("CP");
      for (int j=0;j<jessTimes.size();j++) {
          jessSeries.add(jessTimes.get(j).getFirst(),jessTimes.get(j).getSecond());
          ciaranSeries.add(ciaranTimes.get(j).getFirst(),ciaranTimes.get(j).getSecond());
      }
      
      final XYSeriesCollection dataset = new XYSeriesCollection( );
      dataset.addSeries(jessSeries);
      dataset.addSeries(ciaranSeries);

      JFreeChart xylineChart = ChartFactory.createXYLineChart(
         "", 
         "lambda",
         "Time", 
         dataset,
         PlotOrientation.VERTICAL, 
         true, true, false);
      
      int width = 640;   /* Width of the image */
      int height = 480;  /* Height of the image */ 
      File XYChart = new File( "/home/jessica/Documents/PhD/SubgraphEmbedding/GraphicalResults/"+patternGraph+".jpeg" ); 
      ChartUtilities.saveChartAsJPEG( XYChart, xylineChart, width, height);
   }
   
   //get range of algorithm runtime values (there is one value for each value of
   //lambda) from file and insert the values into a list
    private static List<DoubleIntPair> extractLambdaTimePair(BufferedReader reader) {
        List<DoubleIntPair> lambdaTimePairs = new ArrayList<>();
        DoubleIntPair LambdaTimePair = new DoubleIntPair(-1,-1);
        try {
            String line = reader.readLine();
            
            while (line != null && !line.isEmpty()) {
                
                int i=0;
                String lambdaString = "";
                while (i<line.length()&&Character.isDigit(line.charAt(i))) {
                    lambdaString = lambdaString + line.charAt(i);
                    i++;
                }
                double lambda = Double.parseDouble("2."+lambdaString);
                while (i<line.length()&&!Character.isDigit(line.charAt(i))) {
                    i++;
                }
                String runtimeString = "";
                while (i<line.length()&&Character.isDigit(line.charAt(i))) {
                    runtimeString = runtimeString + line.charAt(i);
                    i++;
                }
                int runtime = Integer.parseInt(runtimeString);
                DoubleIntPair lambdaTimePair = new DoubleIntPair(lambda,runtime);
                lambdaTimePairs.add(lambdaTimePair);
                line = reader.readLine();
            }
            
        } catch (IOException ex) {
            Logger.getLogger(XYLineChart_image.class.getName()).log(Level.SEVERE, null, ex);
        }
        return lambdaTimePairs;
    } 
   
   
   private static BufferedReader getReader(String name) {
        BufferedReader in = null;
        try {
            in = new BufferedReader(new FileReader(new File(name)));
        } catch (FileNotFoundException e) {
            System.out.print("File " + name + " cannot be found");
        }
        return in;
    }
}