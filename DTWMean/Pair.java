/**
 * A class implementing a generic pair of two values.
 *
 * @author Vincent Froese
 * Date: 26 January 2018
 *
 */
public class Pair<T, U> {
        
    public final T t;
    public final U u;

    public Pair(T t, U u) {         
        this.t= t;
        this.u= u;
     }

    public T getFirst() {
	return this.t;
    }

    public U getSecond() {
	return this.u;
    }
}
