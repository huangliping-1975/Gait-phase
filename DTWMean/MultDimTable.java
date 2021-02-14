import java.util.*;

/**
 * A class implementing a multidimensional table as a single array.
 *
 * @author Vincent Froese
 * Date: 9 January 2018
 *
 */
public class MultDimTable<T> {

    private int DIM;
    private int[] SIZE;
    private int[] BASE;
    private int TOTAL_SIZE;
    private ArrayList<T> table;    
    
    public MultDimTable(int[] size) {
	DIM = size.length;
	SIZE = size;
	TOTAL_SIZE = size[0];
	BASE = new int[DIM];
	BASE[0] = 1;
	for (int d=1; d<DIM; d++) {
	    BASE[d] = TOTAL_SIZE;
	    TOTAL_SIZE *= size[d];
	}
	table = new ArrayList<T>(TOTAL_SIZE);
	for (int i=0; i<TOTAL_SIZE; i++) table.add(null);
    }

    private int computeIndex(int[] indices) {
	int i = 0;
	for (int d=0; d<DIM; d++) i += indices[d] * BASE[d];
	return i;
    }
    
    public void put(int[] indices, T data) {
	table.set(computeIndex(indices), data);
    }

    public T get(int[] indices) {
	return table.get(computeIndex(indices));
    }
}
