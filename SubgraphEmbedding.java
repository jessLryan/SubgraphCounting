package subgraphembedding;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Objects;
import java.util.concurrent.TimeUnit;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author jessica
 */
public class SubgraphEmbedding {

    /**
     *
     */
    public static long recCount = 0;
    
    /**
     *
     * @param args
     * @throws java.lang.Exception
     */
    public static void main(String[] args) throws Exception {
        Stats s = new Stats();
        long startTime = System.nanoTime();
        s.startDateTime = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format(new Date());
        s.HFile = args[0];
        s.GFile = args[1];
        BufferedReader HReader = getReader(s.HFile);
        BufferedReader GReader = getReader(s.GFile);
        Graph H = readToHGraph(HReader);
        Graph G = readToGraph(GReader);
        s.GOrder = G.order();
        s.GEdges = G.edgeCount();
        s.GAvgDeg = s.GEdges*2/s.GOrder;
        s.HOrder = H.order();
        s.HEdges = H.edgeCount();
        s.HAvgDeg = s.HEdges*2/s.HOrder;
        Runtime.getRuntime().addShutdownHook(new Thread() {
                @Override
		public void run() {
                    createOutput(s.GFile, s.HFile, s.GOrder, s.GAvgDeg, s.GEdges, s.HOrder, s.HEdges, s.HAvgDeg, s.startDateTime, s.runTime, s.result, s.threshold, s.delta, s.status);   
		}
	    });
        List<Integer> isomorphisms = embed(G, H);
        s.endTime = System.nanoTime();
        s.runTime = TimeUnit.MILLISECONDS.convert(s.endTime - startTime, TimeUnit.NANOSECONDS);
        s.result = Integer.toString(isomorphisms.get(2));
        s.delta = Integer.toString(isomorphisms.get(1));
        s.threshold = Integer.toString(isomorphisms.get(0));
    }
    

    /**
     *
     * @param GFile
     * @param HFile
     * @param GEdges
     * @param GOrder
     * @param GAvgDeg
     * @param HOrder
     * @param HEdges
     * @param HAvgDeg
     * @param startDateTime
     * @param runTime
     * @param result
     * @param threshold
     * @param delta
     * @param status
     */
    public static void createOutput(String GFile, String HFile, int GEdges, int GOrder, int GAvgDeg, int HOrder, int HEdges, int HAvgDeg, String startDateTime, long runTime, String result, String threshold, String delta, String status) {
            System.out.println("Status = "+status);
            System.out.println("Host graph filepath  = "+GFile);
            System.out.println("Host graph order = "+GOrder);
            System.out.println("Host graph edge count = "+GEdges);
            System.out.println("Host graph average degree = "+GAvgDeg);
            System.out.println("Pattern graph filepath  = "+HFile);
            System.out.println("Pattern graph order = "+HOrder);
            System.out.println("Pattern graph edge count = "+HEdges);
            System.out.println("Pattern graph average degree = "+HAvgDeg);
            System.out.println("Start datetime = " + startDateTime);
            System.out.println("runtime = " + runTime);
            //how many of the highest degree vertices we consider to be the
            //"high degree" vertices
            System.out.println("parameter k = " + threshold);
            System.out.println("max degree of remaining vertices = "+delta);
            System.out.println("Number of recursions = " + recCount);
            System.out.println("Class path "+System.getProperty("java.class.path"));
            System.out.println("JRE vendor name "+System.getProperty("java.vendor"));
            System.out.println("JRE version number "+System.getProperty("java.version"));
            System.out.println("Operating system architecture "+System.getProperty("os.arch"));
            System.out.println("Operating system name "+System.getProperty("os.name"));
            System.out.println("Operating system version "+System.getProperty("os.version"));
            System.out.println("User working directory "+System.getProperty("user.dir"));
            System.out.println("User home directory "+System.getProperty("user.home"));
            System.out.println("User account name "+System.getProperty("user.name"));
    }

    //prepares set up for counting copies of H in G. 
    //G is the host graph and H the pattern graph 
    private static List<Integer> embed(Graph G, Graph H) {
        while (1 == 1) {
            if (Thread.currentThread().isInterrupted()) {
                throw new IllegalStateException("time out");
            }
            //add back in if we know H is connected
            //G.removeDegZero();
            List<Integer> results = new ArrayList<>();//contains various result parameters
            int HOrder = H.order();
            double HOrderDouble = (double) HOrder;
            double GOrder = G.order();
            double GMaxDeg = G.maxDeg();
            IntPair kDelta = getThresholdParam(G, GOrder, GMaxDeg, HOrderDouble);
            int chosenK = (int) kDelta.getFirst();
            //delta is the maximum degree of the remaining (i.e. "non high-degree)
            //vetices
            int delta = (int) kDelta.getSecond();
            results.add(chosenK);
            results.add(delta);
            //M contains, for each vertex v in H, the list of vertices in G that v can
            //currently be mapped to
            List<List<Integer>> M = new ArrayList<>();
            for (int l = 0; l < HOrder; l++) {
                List<Integer> mapList = new ArrayList<>();
                M.add(mapList);
            }
           //initialise M to contain, for each v in H, all those vertices in G which 
            //have degree at least the degree of v
            for (int i = 0; i < HOrder; i++) {
                List<Integer> mapSet = M.get(i);
                for (int j = 0; j < GOrder; j++) {
                    if (H.getVertexAtIndex(i).degree() <= G.getVertexAtIndex(j).degree()) {
                        mapSet.add(j);
                    }
                    
                }
                //If any vertex in H cannot be mapped to any vertex in G then there are 
                //no embeddings so return 0
                if (mapSet.isEmpty()) {
                    results.add(0);
                    return results;
                }
            }
            //map high degree vertices first
            int result = 0;
            //get indices of vertices in G with degree above the threshold k
            List<Integer> thresholdVerticesIndices = G.getHighestkIndices(chosenK);
            //mapToThreshold defines the maximum number of vertices in H we try to map to 
            //above threshold vertices
            //either map i=0,...,H.order() vertices of H to the 'above threshold vertices'
            //(in all possible ways) or, if there are fewer above threshold vertices than 
            //vertices in H, then map i=0,...,numberOfThresholdVertices to vertices in H 
            //(in all possible ways) 
            int mapToThreshold = H.order();
            if (H.order() > chosenK) {
                mapToThreshold = chosenK;
            }
            //for all subsets of vertices in H of size i, map these vertices to high
            //degree vertices only, and map the rest of the vertices in H to non-high-degree
            //vertices only
            for (int i = 0; i <= mapToThreshold; i++) {
                //find all sets of i vertices from H
                List<List<Integer>> sets = combination(HOrder, i);
                setH:
                for (List<Integer> set : sets) {
                    List<List<Integer>> MClone = cloneMapSet(M);
                    Graph HClone = H.copy();
                    //update M so that vertices in H with index in set can only
                    //be mapped to high degree vertices in G, and all remaining 
                    //vertices in H can only be mapped to non-high-degree vertices
                    //in G
                    for (int j = 0; j < H.order(); j++) {
                        if (set.contains(j)) {
                            MClone.get(j).clear();
                            MClone.get(j).addAll(thresholdVerticesIndices);
                        } else {
                            MClone.get(j).removeAll(thresholdVerticesIndices);
                            //if any vertex can now not be mapped to any vertex in G, 
                            //then give up on this particular set
                            if (MClone.get(j).isEmpty()) {
                                continue setH;
                            }
                        }
                    }
                    result += mapHighDegree(G, H, MClone, set, 0);
                    
                }
            }
            results.add(result);
            return results;
            
        }
    }
    
    //determine how many of the highest degree vertices in G should be the 
    //"high degree vertices" by optimising theoretical runtime function
    private static IntPair getThresholdParam(Graph G, double GOrder, double GMaxDeg, double HOrderDouble) {
            int delta = -1;
            int oldFunVal = (int) (java.lang.Math.pow(GMaxDeg, HOrderDouble) * java.lang.Math.pow(GOrder, GOrder));
            int chosenK = 0;
            for (int K = 1; K <= java.lang.Math.min(9, GOrder); K++) {
                delta = G.maxDegRemaining(K);
                
                int test = (int) (java.lang.Math.pow(delta, HOrderDouble) * java.lang.Math.pow(K, K));
                double newFunVal = (java.lang.Math.pow(delta, HOrderDouble) * java.lang.Math.pow(K, K));
                if (newFunVal <= 2147483647 && newFunVal < oldFunVal) {
                    chosenK = K;
                    oldFunVal = (int) newFunVal;
                }
            }
            IntPair kDelta = new IntPair(chosenK,delta);
            return kDelta;
    }
    
    
    //mapHighDegree assigns the vertices in H with index in set to high degree vertices
    //in G. It then finds the connected components among the remaining vertices of H and,
    //for each such component, reorders its list L of vertices so that each vertex in L
    //is preceeded by at least one of its neighbours. This is so that we can now call 'map'
    //in the usual way on each component C and it will be FPT since the number of 
    //vertices in G which a given vertex (except the first) in C can be mapped to is bounded
    //by the maximum degree of G
    //G and H are the original host and pattern graph respectively. M contains the list of
    //vertices in G which each vertex in H can be assigned to. set contains the vertices in 
    //H which we are to assign to high degree vertices of G, and k is the number of vertices in H
    //which we have already assigned (to high degree vertices of G) i.e. it keeps track of how many
    //iterations of mapHighDegree we have already performed
    private static double mapHighDegree(Graph G, Graph H, List<List<Integer>> M, List<Integer> set, int k) {
        recCount++;
        double count = 0;
        //if we have mapped all necessary vertices in H to high degree vertices of G, then find the
        //connected components formed by the remaining vertices in H and count the ways in which these
        //components can be embedded in the rest of G using countComponents()
        if (k == set.size()) {
            if (H.order() == set.size()) {
                return 1;
            }
            //a "component" is simply a graph (H) together with a list (M) of
            //vertices (well, their indices) in G to which each vertex in H can
            //be assigned
            Graph HCopy = H.copy();
            Component C = new Component(HCopy, M);
            //remove high degree vertices from C
            Component CNew = C.removeFromSet(set);
            Graph HClone = CNew.getGraph();
            List<List<Integer>> MCopy = CNew.getM();
            //get the connected components of HClone
             List<Component> components = HClone.getComponents(MCopy);
            for (Component D : components) {
                //reorder the vertices in each component so that each vertex
                //is preceeded by at least one of its neighbours
                D.orderVertices();
                //count copies of the component in G
                int compCount = map(G,D.getGraph(),D.getM(),0);
                //if any component has no copies in G, then there can be no
                //copy of H in G with this choice of "high degree set"
                if (compCount == 0) {
                    return 0;
                }
                else if (count == 0) {
                    count = compCount;
                }
                //to get the number of (possibly overlapping) copies of the 
                //components in G, multiply the number of copies of each 
                //individual component in G
                else {
                    count = count*compCount;
                }                
            }
            //for each way to partition the components into "groups", count the 
            //copies of these components in G so that the components in a given
            //group must overlap
            //subtract all of these counts from the original count for the 
            //components
            List<String> list = new ArrayList<>();
            for (int j=0;j<components.size();j++) {
                list.add(Integer.toString(j));
            }
            for (int blocks = 1; blocks <= list.size()-1; ++blocks) {
                for (List<List<String>> partition : new PartitionIterable<>(list, blocks)) {
                    count -= countComponents(G,partition,components);
                }
            }
        } else {
            //otherwise, count the number of ways to assign the next vertex to 
            //each possible high degree vertex of G
            int HIndex = set.get(k);
            Vertex vH = H.getVertexAtIndex(HIndex);
            //indices contains the indices of all vertices in H except those
            //which have already been mapped to high degree vertices. We want
            //to update the mappable sets of these vertices
            List<Integer> indices = new ArrayList<>();
            for (int i = 0; i < M.size(); i++) {
                indices.add(i);
            }
            for (int c = 0; c <= k; c++) {
                indices.remove(set.get(c));
            }
            //assign vH to each possible high degree vertex vG in its mappable set
            assignVG:
            for (int gIndex : M.get(HIndex)) {
                List<List<Integer>> MCopy = cloneMapSet(M);
                //remove vG from all lists in MCopy
                for (Integer j : indices) {
                    //remove vertex in G which has been assigned to vH from 
                    //the mappable sets of vertices in H which have not already
                    //been mapped to high degree vertices
                    MCopy.get(j).remove(Integer.valueOf(gIndex));
                    //if any vertex in H now has an empty mappable set, then we
                    //give up on this particular embedding
                    if (MCopy.get(j).isEmpty()) {
                        continue assignVG;
                    }
                }
                //for each neighbour of vH, remove any vertices in its mappable set which are
                //not neighbours of vG
                for (Vertex hNeigh : vH.getNeighbours()) {
                    int hNeighInd = H.getIndexOfVertex(hNeigh);
                    if (indices.contains(hNeighInd)) {
                        List<Integer> mappable = MCopy.get(hNeighInd);
                        for (Iterator<Integer> itr = mappable.iterator(); itr.hasNext();) {
                            Integer i = itr.next();
                            if (G.getVertexAtIndex(gIndex).getNeighbours().contains(G.getVertexAtIndex(i)) == false) {
                                itr.remove();
                            }
                            //if a neighbour of vH now has an empty mappable set, then we can
                            //give up on this particular embedding
                            if (MCopy.get(hNeighInd).isEmpty()) {
                                continue assignVG;
                            }
                        }
                    }
                }
                count += mapHighDegree(G, H, MCopy, set, k + 1);
            }
        }
        return count;
    }

    //returns a list of all subsets of size k from a set (of integers) of size n
    private static List<List<Integer>> combination(int n, int k) {
        List<List<Integer>> result = new ArrayList<>();
        List<Integer> partialList = new ArrayList<>();
        backtrack(n, k, 0, result, partialList);
        return result;
    }

    
    //used by combination(n,k) to form a subset of size k from n
    private static void backtrack(int n, int k, int startIndex, List<List<Integer>> result, List<Integer> partialList) {
        if (k == partialList.size()) {
            result.add(new ArrayList<>(partialList));
            return;
        }
        for (int i = startIndex; i < n; i++) {
            partialList.add(i);
            backtrack(n, k, i + 1, result, partialList);
            partialList.remove(partialList.size() - 1);
        }
    }

    //used to clone the set M of mappable vertices of H into G

    /**
     *
     * @param M
     * @return
     */
    public static List<List<Integer>> cloneMapSet(List<List<Integer>> M) {
        List<List<Integer>> MCopy = new ArrayList<>();
        for (int i = 0; i < M.size(); i++) {
            List<Integer> mapListCopy = new ArrayList<>();
            MCopy.add(mapListCopy);
            MCopy.get(i).addAll(M.get(i));
        }
        return MCopy;
    }

    
    
    
    //counts copies of remaining components of H in non-high-degree part of G 
    //once vertices in "set" have been mapped to high degree vertices of G
    //"partition" tells us how the components should be "grouped" (the components
    //in a given group must intersect to form a single connected graph). Note that 
    //we may count copies where components in different groups intersect - we 
    //subtract these counts.
    private static double countComponents(Graph G, List<List<String>> partition, List<Component> components) {
        if (components.size()==1) {
            Component C = components.get(0);
            return map(G,C.getGraph(),C.getM(),0);
        }
        //each list in compCollection contains the ways (graphs) in which the
        //components in a group can intersect (one list per group)
        List<List<Component>> compCollection = new ArrayList<>();
        int count = 0;
        for (List<String> part:partition) {
            //if there is only one component in a group, the number of ways in 
            //which the group intersects is just the single component
            if (part.size()==1) {
                List<Component> singleCompList = new ArrayList<>();
                singleCompList.add(components.get((Integer.valueOf(part.get(0)))));
                compCollection.add(singleCompList);
            }
            else {
                List<Component> comps = new ArrayList<>();
                List<Component> compsToIntersect = new ArrayList<>();
                for(String s : part) {
                    compsToIntersect.add(components.get((Integer.valueOf(s))));
                }
                List<Component> interG = intersectingGraphs(compsToIntersect);
                comps.addAll(interG);
                compCollection.add(comps); 
            }
        }
        //if any group of components cannot intersect in any way, there can be
        //no count for them, so return 0
        for (List<Component> compList:compCollection) {
            if (compList.isEmpty()) {
                return 0;
            }
        }        
        //for each way a group of components can intersect, we must count copies of
        //this intersection in G together with each way that all other groups can intersect
        //and be mapped to G
        //to do this, we must obtain all combinations of intersections in groups (one 
        //intersection from each group at a time) using listPermutations
        int[] limits = new int[compCollection.size()];
        for (int i=0;i<compCollection.size();i++) {
            limits[i] = compCollection.get(i).size()-1;
        }
        int[] indices = new int[compCollection.size()];
        List<int[]> indexList = listPermutations(indices,limits, false);
        //for each collection of intersecting components (with one from
        //each group in each case), count copies of these intersections in G
        //using getCOunt (which also removes all cases where the intersections
        //also intersect each other in G)
        for (int[] index :indexList) {
            List<Component> compList = new ArrayList<>();
            for (int i=0;i<compCollection.size();i++) {
                compList.add(compCollection.get(i).get(index[i]));
            }
            int getCount = getCount(G,compList);
            count+= getCount;
        }
        return count;
    }
    
    
    //returns all sequences of limits.size() numbers where position i in the 
    //sequence can have value at most limits[i] (and starts from value 0)
    private static List<int[]> listPermutations(int[] indices, int[] limits, boolean firstIsArbitrary) {
        //addTo contains all the sequences
        List<int[]> addTo = new ArrayList<>();
        int[] indicesCopy = new int[indices.length];
        int lowerBound = 0;
        if (firstIsArbitrary) {
            indicesCopy[0]=0;
            lowerBound = 1;
        }
        System.arraycopy(indices, 0, indicesCopy, 0, indices.length);
        //add current sequence
        addTo.add(indices);
        int size = limits.length;
        boolean foundIndex = false;
        int checkIndex = size-1;
        //works by increasing the rightmost feasible digit in the sequence 
        //(i.e. the rightmost digit which is below the corresponding value in
        //limits). If a digit is increased which is not the final digit in the 
        //sequence, all digits to the right are set to zero before proceeding
        while (!foundIndex) {
            if (limits[checkIndex]==indicesCopy[checkIndex]) {
                if (checkIndex==lowerBound) {
                    //if the index we are at has limit zero and is the leftmost
                    //index, we  are done
                    foundIndex = true;
                }
                else {
                    checkIndex--;
                }
            }
            //if we can increase the digit at this index, do so, and reset all
            //the digits further right in the sequence to 0
            else if (indicesCopy[checkIndex]<limits[checkIndex]) {
                indicesCopy[checkIndex]++;
                foundIndex = true;
                for (int j=checkIndex+1;j<size-1;j++) {
                    indicesCopy[j]=0;
                }
                //this is the line I temporarily delete
                addTo.addAll(listPermutations(indicesCopy,limits, firstIsArbitrary));
            }
            else {
                //if we can't further increase the digit at this index, move 
                //further to the left along the sequence
                if (checkIndex==lowerBound) {
                    //if we are at the leftmost digit in the sequence, we are 
                    //done
                    foundIndex = true;
                }
                checkIndex--;
            }
        }
        return addTo;       
    }
        
 
    //count copies of components in G and subtract intersections
    private static int getCount(Graph G, List<Component> components) {
        int count = 1;
        for (Component C:components) {
            C.orderVertices();
            int compCount = map(G,C.getGraph(),C.getM(),0);
            if (compCount == 0) {
                return 0;
            }
            else if (count == 0) {
                count = compCount;
            }
            else {
                count = count*compCount;
            }
        }
        //find ways components can intersect and subtract count
        List<String> list = new ArrayList<>();
        for (int j=0;j<components.size();j++) {
            list.add(Integer.toString(j));
        }
        if (components.size()>1) {
            for (int blocks = 1; blocks <= list.size()-1; ++blocks) {
                for (List<List<String>> partition : new PartitionIterable<>(list, blocks)) {
                    boolean somePartGreaterThanOne = false;
                    for (List<String> part:partition) {
                        if (part.size()>1) {
                            somePartGreaterThanOne = true;
                        }
                    }
                    if (somePartGreaterThanOne) {
                            count -= countComponents(G,partition,components);
                    }
                }
            }
        }
        return count;
    }

    
    //returns ways in which compsToIntersect can intersect to form a connected graph
    //this is done by creating a graph G with number of vertices equal to |C_0|
    //plus the sum of C_i-1 for each component C_i in comps to intersect (except C_0),
    //and then counting the ways in which the components can be mapped in G (taking repetitions into account)
    private static List<Component> intersectingGraphs(List<Component> compsToIntersect) {
        int GOrder = compsToIntersect.get(0).getGraph().order();
            for (int i=1;i<compsToIntersect.size();i++) {
                GOrder += compsToIntersect.get(i).getGraph().order()-1;
            }
        Graph G = new Graph(GOrder);
        List<List<List<Integer>>> mappedlists = new ArrayList<>();
        //create mapList for each component (each vertex can initially be 
        //mapped to any vertex in the created graph) and find copies of the
        //component in G
        List<Integer> firstCompCurrentList = new ArrayList<>();
        List<List<Integer>> firstCompTotalList = new ArrayList<>();
        
        int firstCompSize = compsToIntersect.get(0).getGraph().order();
        for (int k=0;k<firstCompSize;k++) {
            firstCompCurrentList.add(k);
        }
        firstCompTotalList.add(firstCompCurrentList);
        mappedlists.add(firstCompTotalList);
        for (int j=1;j<compsToIntersect.size();j++) {
            Component currentComp = compsToIntersect.get(j);
            Graph currentH = currentComp.getGraph();
            List<List<Integer>> currentM = new ArrayList<>();
            List<Integer> currentMList = new ArrayList<>();
            for (int i = 0; i < G.order(); i++) {
                currentMList.add(i);
            }
           for (int i = 0; i < currentH.order(); i++) {
               List<Integer> newList = new ArrayList<>();
               newList.addAll(currentMList);
               currentM.add(i, newList);
           }
           //find the ways in which currentComp can be mapped to the vertices in 
           //G
           List<Integer> currentList = new ArrayList<>();
           List<List<Integer>> totalList = new ArrayList<>();
           returnVertices(G, currentH, currentM, 0, currentList, totalList);
           mappedlists.add(totalList);
        }
        //reorder components and corresponding mapped lists so that the
        //biggest component comes first. This is so that when we try to place
        //the components in the "constructed graph", we can assign the biggest
        //component only once (since the position of the first component is
        //arbitrary, only the relative positions of the other components matter)
        int biggestCompIndex = -1;
        int size = 0;
        for (Component C:compsToIntersect) {
            if (C.getGraph().order()>size) {
                biggestCompIndex = compsToIntersect.indexOf(C);
            }
        }
        Component bigComp = compsToIntersect.get(biggestCompIndex);
        Component zeroComp = compsToIntersect.get(0);
        compsToIntersect.set(0, bigComp);
        compsToIntersect.set(biggestCompIndex, zeroComp);
        List<List<Integer>> bigList = mappedlists.get(biggestCompIndex);
        List<List<Integer>> zeroList = mappedlists.get(0);
        mappedlists.set(0, bigList);
        mappedlists.set(biggestCompIndex, zeroList);
        //need to consider all possible combinations of ways each component can
        //be mapped into G (use listPermutations)
        int[] limits = new int[mappedlists.size()];
        int[] indices = new int[mappedlists.size()];
        for (int i=0;i<limits.length;i++) {
            limits[i] = mappedlists.get(i).size()-1;
        }
        List<int[]> indicesList = listPermutations(indices,limits, false); 
        List<int[]> myList = new ArrayList<>();
        int index = 0;
        //get and return intersection graphs using obtainIntersectionGraphs and the sequence lists 
        //created using listPermutations. "index" tells us which list generated
        //using listPermutations we use in a given iteration of obtainIntersectionGraphs
        return obtainIntersectGraphs(compsToIntersect,mappedlists,indicesList,index,myList);
    }

    
    //returns the ways in which a collection of graphs can intersect to form a 
    //connected graph after "mappedlists" has been created in the 
    //intersectingGraphs method
    private static List<Component> obtainIntersectGraphs(List<Component> compsToAdd, List<List<List<Integer>>> mappedlists, List<int[]> indices, int index, List<int[]> dupCheck) {
        List<Component> returnComps = new ArrayList<>();
        if (compsToAdd.size()==1) {
            returnComps.add(compsToAdd.get(0));
            return returnComps;
        }
        List<Integer> allMapsets = new ArrayList<>();
        //get the specific list of vertices in G (the original constructed graph)
        //to which the vertices of each component are mapped to create this 
        //intersecting graph
        int[] indexList = indices.get(index);
        for (int i=0;i<indexList.length;i++) {
            List<Integer> mapSet = new ArrayList<>();
            mapSet.addAll(mappedlists.get(i).get(indexList[i]));
            allMapsets.addAll(mapSet);
        }
        //each component must intersect with at least one other component. To
        //check for this, we check, for each component, that at least one
        //vertex is mapped to is also mapped to by a vertex from a different 
        //component
        boolean verticesOverlap = true;
        int firstIndex = 0;
        comps: for (Component C:compsToAdd) {
            int endIndex = firstIndex+C.getGraph().order();
            List<Integer> sublist = allMapsets.subList(firstIndex, endIndex);
            boolean someVertexOverlaps = true;
            for (Integer s:sublist) {
                for (int k=0;k<allMapsets.size();k++) {
                    if (Objects.equals(allMapsets.get(k), s)&&(k<firstIndex||k>endIndex-1)) {
                        firstIndex+=C.getGraph().order();
                        continue comps;
                    }
                }
            }
            verticesOverlap = false;
            break;
        }
        if (verticesOverlap) {
            //create intersecting graph
            Graph constructedGraph = new Graph();
            List<List<Integer>> MCloned = new ArrayList<>();
            for (Component C : compsToAdd) {
                Graph F = C.getGraph();
                //initially graph is just made up of original components all
                //stuck together, then vertices are "merged"
                constructedGraph.addNewVertices(F.getVertices());
                for (List<Integer> list:C.getM()) {
                    List<Integer> newList = new ArrayList<>();
                    newList.addAll(list);
                    MCloned.add(newList);
                }
            }
            //dupCheck stores lists of where all vertices were mapped to the 
            //constructed graph in each possibility. We then compare this mapping 
            //to each mapping stored in dupCheck to see if we have already 
            //created an identical intersection graph from a different mapping
            int graphOrder = constructedGraph.order();
            int[] dupCheckArr = new int[graphOrder];
            for (int i=0;i<graphOrder;i++) {
                //dupCheckArr initially "assumes" all vertices are mapped to 
                //distinct vertices, then updates later
                dupCheckArr[i] = i;
            }
            //if two vertices are mapped to the same vertex they are given the
            //same value in dupCheckArr
            for (int b=0;b<allMapsets.size();b++) {
                Integer mapped = allMapsets.get(b);
                for (int c=b+1;c<allMapsets.size();c++) {
                    if (Objects.equals(mapped, allMapsets.get(c))) {
                        dupCheckArr[c]=b;
                    }
                }
            }
            //if we have already created an identical intersection graph from a
            //different mapping then we discard this mapping
            boolean use = true;
            for (int[] arr:dupCheck) {
                if (Arrays.equals(dupCheckArr, arr)) {
                    use=false;
                }
            }
            //otherwise, we begin to construct the intersection graph
            if (use) {
                dupCheck.add(dupCheckArr);
                //stores the maplists that we will "remove" due to merged vertices 
                List<Integer> deleteFromMClone = new ArrayList<>();
                //stores which vertices are mapped to the same vertex (used to
                //update MCloned)
                List<List<Integer>> sameAs = new ArrayList<>();
                for (int i=0;i<graphOrder;i++) {
                    List<Integer> sameAsList = new ArrayList<>();
                    sameAs.add(sameAsList);
                }
                //if multiple vertices are mapped to a vertex then merge them
                for (int o=0;o<graphOrder;o++) {
                    boolean found = false;
                    int indexx = 0;
                    while(!found&&indexx<allMapsets.size()) {
                        if (o==allMapsets.get(indexx)) {
                            found = true;
                            //look for other vertices that are mapped to vertex o
                            for (int otherIndex=indexx+1;otherIndex<allMapsets.size();otherIndex++) {
                                if (o==allMapsets.get(otherIndex)) {
                                    deleteFromMClone.add(otherIndex);
                                    sameAs.get(indexx).add(otherIndex);
                                    //delete all neighbours of vertex at indexx
                                    //that are not meighbours of vertex at 
                                    //otherIndex
                                    constructedGraph.mergeVertices(indexx, otherIndex);
                                    constructedGraph.getVertexAtIndex(otherIndex).setDeleteFlag(true);
                                    //vertex at otherIndex is going to be deleted,
                                    //so delete any adjacencies
                                    constructedGraph.getVertexAtIndex(otherIndex).deleteNeighbours();
                                    
                                    
                                }                            
                            }
                        }
                        else{
                            indexx++;
                        }
                    }                
                }
                //create new updated version of maplist MCloned after merging 
                //vertices
                List<List<Integer>> MClonedCloned = new ArrayList<>();
                for (int i=0;i<graphOrder;i++) {
                    List<Integer> newList = new ArrayList<>();
                    for (Integer j :MCloned.get(i)) {
                        boolean include = true;
                        for (int k:sameAs.get(i)) {
                            if (!MCloned.get(k).contains(j)) {
                                include = false;
                            }
                        }
                        if (include) {
                            newList.add(j);
                        }
                    }
                    MClonedCloned.add(newList);
                }
                //only include non "duplicate" vertices from MClonedMCloned
                List<List<Integer>> MClonedClonedCloned = new ArrayList<>();
                for (int i=0;i<MClonedCloned.size();i++) {
                    if (!deleteFromMClone.contains(i)) {
                        MClonedClonedCloned.add(MClonedCloned.get(i));
                    }
                }
                //delete duplicate vertices from constructedGraph
                Iterator<Vertex> iter = constructedGraph.getVertices().iterator();
                while (iter.hasNext()) {
                    Vertex v = iter.next();
                    if (v.getDeleteFlag()) {
                        iter.remove();
                    }
                }
                //if any list in MClonedClonedCloned is empty, there will be no
                //copies of this intersection graph in G, so don't include it
                boolean add = true;
                for (List<Integer> MList :MClonedClonedCloned) {
                    if (MList.isEmpty()) {
                        add = false;
                    }
                }
                if (add) {
                    Component constructedComp = new Component(constructedGraph, MClonedClonedCloned);
                    Graph consG = constructedComp.getGraph();
                    returnComps.add(constructedComp);
                }
            }
            
        }
        if (index<indices.size()-1) {
            returnComps.addAll(obtainIntersectGraphs(compsToAdd, mappedlists, indices,index+1,dupCheck));
        }
        return returnComps;
    }
    
    
    
    
    //create a deep copy of a list of integers

    /**
     *
     * @param intList
     * @return
     */
    public static List<Integer> copyListInt(List<Integer> intList) {
        List<Integer> newList = new ArrayList<>();
        intList.stream().forEach((i) -> {
            newList.add(i);
        });
        return newList;
    }

    //find each copy of H in G (using the same idea as the method map()) and, 
    //for each copy, return the vertices of G to which each vertex of H is 
    //mapped. Each list is formed using currentList. All the lists are contained 
    //in totalList
    private static void returnVertices(Graph G, Graph H, List<List<Integer>> M, int HIndex, List<Integer> currentList, List<List<Integer>> totalList) {
        recCount++;
        Vertex vH = H.getVertexAtIndex(HIndex);
        //assign vertex in G with index gIndex to vertex in H with index HIndex
        assignVG:
        for (int gIndex : M.get(HIndex)) {
            List<List<Integer>> MCopy = cloneMapSet(M);
            List<Integer> currentListCopy = copyListInt(currentList);
            currentListCopy.add(gIndex);
            if (HIndex < H.order() - 1) {
                for (int j = HIndex + 1; j < M.size(); j++) {
                    //remove vertex in G which has been assigned to vH from all lists in MCopy
                    MCopy.get(j).remove(Integer.valueOf(gIndex));
                    if (MCopy.get(j).isEmpty()) {
                        continue assignVG;
                    }
                }
            }
            //for each neighbour of vH, remove any vertices in its mappable set which are
            //not neighbours of vG
            for (Vertex hNeigh : vH.getNeighbours()) {
                int hNeighInd = H.getIndexOfVertex(hNeigh);
                //only need to consider mappable set of neighbours with index at least 
                //that of vH since all other neighbours have already been assigned to
                //some vertex in G
                if (hNeighInd > HIndex) {
                    List<Integer> mappable = MCopy.get(hNeighInd);
                    for (Iterator<Integer> itr = mappable.iterator(); itr.hasNext();) {
                        Integer i = itr.next();
                        if (G.getVertexAtIndex(gIndex).getNeighbours().contains(G.getVertexAtIndex(i)) == false) {
                            itr.remove();
                        }
                        //if a neighbour of vH now has an empty mappable set, then we can
                        //give up on this particular embedding
                        if (MCopy.get(hNeighInd).isEmpty()) {
                            continue assignVG;
                        }
                    }
                }
            }
            //if there are vertices in H which have not yet been assigned, recurse map().
            //Else, we have found one embedding
            if (HIndex < H.order() - 1) {
                returnVertices(G, H, MCopy, HIndex + 1, currentListCopy, totalList);
            } else {
                totalList.add(currentListCopy);
            }
        }
    }
       
    
    
    //in each recursion, map() assigns a vertex vH in H to each possible vertex in  
    //its mappable set in M (and updates M_j for each neighbour v_j of vH)
    private static int map(Graph G, Graph H, List<List<Integer>> M, int HIndex) {
        recCount++;
        int embeddings = 0;
        Vertex vH = H.getVertexAtIndex(HIndex);
        //assign vertex in G with index gIndex to vertex in H with index HIndex
        assignVG:
        for (int gIndex : M.get(HIndex)) {
            List<List<Integer>> MCopy = cloneMapSet(M);
            for (int j = HIndex + 1; j < M.size(); j++) {
                //remove vertex in G which has been assigned to vH from all lists in MCopy
                MCopy.get(j).remove(Integer.valueOf(gIndex));
                if (MCopy.get(j).isEmpty()) {
                    continue assignVG;
                }
            }
            //for each neighbour of vH, remove any vertices in its mappable set which are
            //not neighbours of vG
            for (Vertex hNeigh : vH.getNeighbours()) {
                int hNeighInd = H.getIndexOfVertex(hNeigh);
                //only need to consider mappable set of neighbours with index at least 
                //that of vH since all other neighbours have already been assigned to
                //some vertex in G
                if (hNeighInd > HIndex) {
                    List<Integer> mappable = MCopy.get(hNeighInd);
                    for (Iterator<Integer> itr = mappable.iterator(); itr.hasNext();) {
                        Integer i = itr.next();
                        if (G.getVertexAtIndex(gIndex).getNeighbours().contains(G.getVertexAtIndex(i)) == false) {
                            itr.remove();
                        }
                        //if a neighbour of vH now has an empty mappable set, then we can
                        //give up on this particular embedding
                        if (MCopy.get(hNeighInd).isEmpty()) {
                            continue assignVG;
                        }
                    }
                }
            }
            //if there are vertices in H which have not yet been assigned, recurse map().
            //Else, we have found one embedding
            if (HIndex < H.order() - 1) {
                embeddings += map(G, H, MCopy, HIndex + 1);
            } else {
                embeddings++;
            }
        }
        return embeddings;
    }

    public static BufferedReader getReader(String name) {
        BufferedReader in = null;
        try {
            in = new BufferedReader(new FileReader(new File(name)));
        } catch (FileNotFoundException e) {
            System.out.print("File " + name + " cannot be found");
        }
        return in;
    }

    //converts Ciaran's graph format into a graph
    private static Graph readToGraph(BufferedReader reader) {
        Graph G = new Graph();
        try {
            //the first line of the file lists the number of vertices in the graph
            String line = reader.readLine();
            for (int k = 0; k < Integer.parseInt(line); k++) {
                    Vertex u = new Vertex();
                    G.addVertex(u);
            }
            line = reader.readLine();
            int vertexIndex = -1;
            while (line != null&&!line.isEmpty()) {
                if (line.indexOf("newVertex")>=0) {
                    vertexIndex++;
                    Vertex v = G.getVertexAtIndex(vertexIndex);
                    int i = 10;
                    if (line.length() > 10) {
                        while (Character.isDigit(line.charAt(i))) {
                            i++;
                        }
                        for (int j = i; j < line.length(); j++) {
                            String ch = "";
                            ch += line.charAt(j);
                            if (Character.isDigit(line.charAt(j))) {
                                while (j + 1 < line.length() && Character.isDigit(line.charAt(j + 1))) {
                                    ch += line.charAt(j + 1);
                                    j++;
                                }
                                int dig = Integer.parseInt(ch);
                                Vertex u = G.getVertexAtIndex(dig);
                                if (vertexIndex != dig && !u.isAdjacent(v)) {
                                    v.makeAdjacent(u);
                                    u.makeAdjacent(v);
                                }
                            }
                        }
                    } 
                                        
                }
                else {
                    Vertex v = G.getVertexAtIndex(vertexIndex);
                    int i = 0;
                    if (line.length() > 1) {
                        for (int j = i; j < line.length(); j++) {
                            String ch = "";
                            ch += line.charAt(j);
                            if (Character.isDigit(line.charAt(j))) {
                                while (j + 1 < line.length() && Character.isDigit(line.charAt(j + 1))) {
                                    ch += line.charAt(j + 1);
                                    j++;
                                }
                                int dig = Integer.parseInt(ch);
                                Vertex u = G.getVertexAtIndex(dig);
                                if (vertexIndex != dig && !u.isAdjacent(v)) {
                                    v.makeAdjacent(u);
                                    u.makeAdjacent(v);
                                }
                            }
                        }
                    } 
                }
                line = reader.readLine();
            }
        } catch (IOException e) {
            System.out.print("I/O Error");
            System.exit(0);
        }
        return G;
    }

    //converts Ciaran's graph format into a graph
    private static Graph readToHGraph(BufferedReader reader) {
        Graph G = new Graph();
        try {
            //the first line of the file lists the number of vertices in the graph
            String line = reader.readLine();
            for (int k = 0; k < Integer.parseInt(line); k++) {
                    Vertex u = new Vertex();
                    G.addVertex(u);
            }
            line = reader.readLine();
            int vertexIndex = 0;
            while (line != null&&!line.isEmpty()) {
                Vertex v = G.getVertexAtIndex(vertexIndex);
                    
                    int i = 0;
                    if (line.length() > 1) {
                        while (Character.isDigit(line.charAt(i))) {
                            i++;
                        }
                        for (int j = i; j < line.length(); j++) {
                            String ch = "";
                            ch += line.charAt(j);
                            if (Character.isDigit(line.charAt(j))) {
                                while (j + 1 < line.length() && Character.isDigit(line.charAt(j + 1))) {
                                    ch += line.charAt(j + 1);
                                    j++;
                                }
                                int dig = Integer.parseInt(ch);
                                Vertex u = G.getVertexAtIndex(dig);
                                if (vertexIndex != dig && !u.isAdjacent(v)) {
                                    v.makeAdjacent(u);
                                    u.makeAdjacent(v);
                                }
                            }
                        }
                    } 
                line = reader.readLine();
                vertexIndex++;
                
            }
        } catch (IOException e) {
            System.out.print("I/O Error");
            System.exit(0);
        }
        return G;
    }

    //gets order and of first index of graph from Jess Enright graph file 
    private static int getOrder(BufferedReader reader) {
        int order = 0;
        try {
            String line = reader.readLine();
            while (line != null) {
                int i = 0;
                String vertex1 = "";
                while (i < line.length() && line.charAt(i) != ',' && line.charAt(i) != ' ') {
                    vertex1 += line.charAt(i);
                    i++;
                }
                if (Integer.parseInt(vertex1) > order) {
                    order = Integer.parseInt(vertex1);
                }
                i++;
                String vertex2 = "";
                while (i < line.length()) {
                    vertex2 += line.charAt(i);
                    i++;
                }
                if (Integer.parseInt(vertex2) > order) {
                    order = Integer.parseInt(vertex2);
                }
                line = reader.readLine();
            }
        } catch (IOException e) {
            System.out.print("I/O Error");
            System.exit(0);
        }
        return order;
    }

}