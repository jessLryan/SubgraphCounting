/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package subgraphembedding;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
/**
 *
 * @author jessica
 * @param <T>
 */
public class PartitionIterable<T> implements Iterable<List<List<T>>> {

    private final List<T> elements = new ArrayList<>();
    private final int blocks;

    /**
     *
     * @param elements
     * @param blocks
     */
    public PartitionIterable(List<T> elements, int blocks) {
        checkNumberOfBlocks(blocks, elements.size());
        this.elements.addAll(elements);
        this.blocks = blocks;
    }

    @Override
    public Iterator<List<List<T>>> iterator() {
        return new PartitionIterator<>(elements, blocks);
    }

    private void checkNumberOfBlocks(int blocks, int numberOfElements) {
        if (blocks < 1) {
            throw new IllegalArgumentException(
                    "The number of blocks should be at least 1, received: " +
                    blocks);
        }

        if (blocks > numberOfElements) {
            throw new IllegalArgumentException(
                    "The number of blocks should be at most " +
                    numberOfElements + ", received: " + blocks);
        }   
    }

    private static final class PartitionIterator<T> 
    implements Iterator<List<List<T>>> {

        private List<List<T>> nextPartition;
        private final List<T> elements = new ArrayList<>();
        private final int blocks;

        private final int[] s;
        private final int[] m;
        private final int n;

        PartitionIterator(List<T> elements, int blocks) {
            this.elements.addAll(elements);
            this.blocks = blocks;
            this.n = elements.size();
            s = new int[n];
            m = new int[n];
                //for (int i = 0; i < n - blocks + 1; ++i) {
                //    s[i] = 0;
                //    m[i] = 0;
                //}
                for (int i = n - blocks + 1; i < n; ++i) {
                    s[i] = m[i] = i - n + blocks;
                }
                loadPartition();
        }

        @Override
        public boolean hasNext() {
            return nextPartition != null;
        }   

        @Override
        public List<List<T>> next() {
            if (nextPartition == null) {
                throw new NoSuchElementException("No more partitions left.");
            }

            List<List<T>> partition = nextPartition;
            generateNextPartition();
            return partition;
        }

        private void loadPartition() {
            nextPartition = new ArrayList<>(blocks);

            for (int i = 0; i < blocks; ++i) {
                nextPartition.add(new ArrayList<>());
            }

            for (int i = 0; i < n; ++i) {
                nextPartition.get(s[i]).add(elements.get(i));
            }
        }

        private void generateNextPartition() {
            for (int i = n - 1; i > 0; --i) {
                if (s[i] < blocks - 1 && s[i] <= m[i - 1]) {
                    s[i]++;
                    m[i] = Math.max(m[i], s[i]);

                    for (int j = i + 1; j < n - blocks + m[i] + 1; ++j) {
                        s[j] = 0;
                        m[j] = m[i];
                    }

                    for (int j = n - blocks + m[i] + 1; j < n; ++j) {
                        s[j] = m[j] = blocks - n + j;
                    }

                    loadPartition();
                    return;
                }
            }

            nextPartition = null;
        }
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {
        List<String> list = Arrays.asList("1", "2");
        int row = 1;
    
        for (int blocks = 1; blocks <= list.size(); ++blocks) {
            for (List<List<String>> partition : new PartitionIterable<>(list, blocks)) {
                System.out.printf("%2d: %s\n", row++, partition);
    }
        }
    }
}
