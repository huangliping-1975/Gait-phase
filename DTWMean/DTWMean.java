import java.util.*;

/**
 * A class implementing several methods to compute dtw-means for a given sample of time series.
 *
 * @author Vincent Froese
 * Date: 19 February 2018
 *
 */
public class DTWMean {

    /**
     * The main method runs the algorithms and prints out all solutions.
     */
    public static void main(String[] args) {

        //input time series
	double[][] X = {{1.0, 10.0, 0.0, 0.0, 4.0},
			{0.0, 2.0, 10.0, 0.0, 0.0},
			{0.0, 0.0, 0.0, 10.0, 0.0}};
	
	double[] weights = {1.0/3.0, 1.0/3.0, 1.0/3.0};
	
        Pair<Double, ArrayList<ArrayList<Double>>> result = DTWMean(X, weights);
	System.out.println(result.getSecond().size() + " condensed weighted mean(s) with Frechet value F = " + result.getFirst());
	for (ArrayList<Double> mean : result.getSecond()) {
	    System.out.println(mean);
	}	
	
	int L=10;
	Pair<Double[], ArrayList<ArrayList<ArrayList<Double>>>> all = fixedLengthDTWMean(X, weights, L);
	for (int l=0; l<L; l++) {
	    System.out.println(all.getSecond().get(l).size() + " length-" + (l+1) + " mean(s) with Frechet value F = " + all.getFirst()[l]);
	    for (ArrayList<Double> mean : all.getSecond().get(l)) {
		System.out.println(mean);
	    }
	}
    }
    
    
    /**
     * Method computing a list of all condensed weighted dtw-means of the input time series X.
     * Returns a pair containing the F-value and a list of means.
     *
     * Parameter w: Array of sample weights between 0.0 and 1.0 summing up to 1.0.
     */
    public static Pair<Double, ArrayList<ArrayList<Double>>> DTWMean(double[][] X, double[] w) {
	int k = X.length;
	int[] N = new int[k];
	
        int[] I = new int[k];
        int[] J = new int[k];
	int[] U = new int[k];
	for (int i=0; i<k; i++) {
	    N[i] = X[i].length;
	}

	//initialize DP tables
	MultDimTable<Double> C = new MultDimTable<Double>(N);       //table storing optimal F-values
	
        MultDimTable<ArrayList<ArrayList<Double>>> Z = new MultDimTable<ArrayList<ArrayList<Double>>>(N);                        //table storing lists of means
	
	double[] res = computeMuSigma(X, w, I, I);
	double mu = res[0];
	double sigma = res[1];
	ArrayList<Double> z = new ArrayList<Double>();
	z.add(mu);
	ArrayList<ArrayList<Double>> meanList = new ArrayList<ArrayList<Double>>();
	meanList.add(z);
	C.put(I, sigma);
	Z.put(I, meanList);

	//fill tables by iterating over all I indices i_1 = 0,...,n_1-1, i_2 =...
	while (true) {
	    //increment I indices
	    boolean carry = true;
	    int d = 0;
	    while (carry && d < k) {
		if (I[d] == N[d] - 1) {
		    I[d] = 0;
		    d++;
		} else {
		    I[d]++; 
		    carry = false;
		}
	    }
	    if (carry) break; //finished I iteration

	    J = new int[k];	    	    

	    res = computeMuSigma(X, w, I, J);
	    mu = res[0];
	    sigma = res[1];
	    C.put(I, sigma);
	    z = new ArrayList<Double>();
	    z.add(mu);
	    meanList = new ArrayList<ArrayList<Double>>();
	    meanList.add(z);
	    Z.put(I, meanList);
		    
	    //iterate over all J indices j_1 = 0,...,i_1, j_2 = ...
	    while (true) {
	        //increment J indices
	        carry = true;
	        d = 0;
		while (carry && d < k) {
		    if (J[d] == I[d]) {
			J[d] = 0;
			d++;
		    } else {
			J[d]++; 
			carry = false;
		    }
		}
		if (carry) break; //finished J iteration
		
		//recursively compute optimal mean
		res = computeMuSigma(X, w, I, J);
		mu = res[0];
		sigma = res[1];
		
		//initialize U indices
		for (int i=0; i<k; i++) {
		    if (J[i] == 0) { // J index already zero
			U[i] = 0;
		    } else { 
			U[i] = J[i]-1;
		    }
		}
		    
		double c = Double.POSITIVE_INFINITY;
		    
		//iterate over all indices u_i in {j_i-1,j_i}
		while (true) {
		    double t = C.get(U);
		    if (t < c) {
			c = t;
			meanList = new ArrayList<ArrayList<Double>>();
			for (ArrayList<Double> x : Z.get(U))
			    meanList.add(new ArrayList<Double>(x));
		    } else if (t == c) {
			for (ArrayList<Double> x : Z.get(U)) {
			    if (!meanList.contains(x))
				meanList.add(new ArrayList<Double>(x));
			}
		    }
			
		    // increment U indices
		    carry = true;
		    d = 0;
		    while (carry && (d < k)) {
			if (U[d] == J[d]) {
			    if (U[d] > 0)
				U[d]--;
			    d++;
			} else {
			    U[d]++; 
			    carry = false;
			}
		    }
		    if (carry) break; //finished U iteration
		}

		sigma += c;
		if (sigma < C.get(I)) {
		    C.put(I, sigma);
		    ArrayList<ArrayList<Double>> Q = new ArrayList<ArrayList<Double>>();
		    for (ArrayList<Double> x : meanList) {
			if (mu != x.get(x.size()-1)) //only keep condensed means
			    x.add(mu);
			if (!Q.contains(x))
			    Q.add(x);
		    }
		    Z.put(I, Q);
		} else if (sigma == C.get(I)) {
		    for (ArrayList<Double> x : meanList) {
			if (mu != x.get(x.size()-1)) //only keep condensed means
			    x.add(mu);
			if (!Z.get(I).contains(x))
			    Z.get(I).add(new ArrayList<Double>(x));
		    }
		}
	    }
	}
	int[] r = new int[k];
	for (int i=0; i<k; i++) r[i] = N[i]-1;
	return new Pair<Double, ArrayList<ArrayList<Double>>>(C.get(r), Z.get(r));
    }
    

    /**
     * Method computing a list containing all weighted dtw-means of the input time series X for all fixed lengths from 1 to maxLength.
     * Returns a pair containing the F-values and an array of mean lists.
     *
     * Parameter w: Array of sample weights between 0.0 and 1.0 summing up to 1.0.
     * Parameter maxLength: Maximum length of a mean.
     */
    public static Pair<Double[], ArrayList<ArrayList<ArrayList<Double>>>> fixedLengthDTWMean(double[][] X, double[] w, int maxLength) {
	int k = X.length;
	int[] N = new int[k+1];
	for (int i=0; i<k; i++) {
	    N[i] = X[i].length;
	}
	N[k] = maxLength;

	int[] I = new int[k];
	int[] J = new int[k];
	int[] U = new int[k];
	int[] H = new int[k+1];

	//initialize DP tables
	
	//table storing optimal F-values
	MultDimTable<Double> C = new MultDimTable<Double>(N);
	//table storing lists of means
        MultDimTable<ArrayList<ArrayList<Double>>> Z = new MultDimTable<ArrayList<ArrayList<Double>>>(N);
	
	//fill tables
	//iterate over all lengths L = 1,...,maxLength
	for (int L=1; L<=maxLength; L++) {

	    H[k] = L-1;
	    for (int i=0; i<k; i++) {
		I[i] = 0;
		H[i] = 0;
	    }
	    double[] res = computeMuSigma(X, w, I, I);
	    double mu = res[0];
	    double sigma = res[1];
	    ArrayList<Double> z = new ArrayList<Double>();
	    for (int a=0; a<L; a++) z.add(mu);
	    ArrayList<ArrayList<Double>> meanList = new ArrayList<ArrayList<Double>>();
	    meanList.add(z);
	    C.put(H, sigma*L);
	    Z.put(H, meanList);
	    
	    //iterate over all I indices i_1 = 0,...,n_1-1, i_2 =...
	    while (true) {
		//increment I indices
		boolean carry = true;
		int d = 0;
		while (carry && d < k) {
		    if (I[d] == N[d] - 1) {
			I[d] = 0;
			H[d] = 0;
			d++;
		    } else {
			I[d]++;
			H[d]++;
			carry = false;
		    }
		}
		if (carry) break; //finished I iteration

		J = new int[k];
		
		if (L==1) {
		    res = computeMuSigma(X, w, I, J);
		    mu = res[0];
		    sigma = res[1];
		    C.put(H, sigma);
		    z = new ArrayList<Double>();
		    z.add(mu);
		    meanList = new ArrayList<ArrayList<Double>>();
		    meanList.add(z);
		    Z.put(H, meanList);
		} else {
		    C.put(H, Double.POSITIVE_INFINITY);
		
		    //iterate over all J indices j_1 = 0,...,i_1, j_2 = ...
		    while (true) {

			//recursively compute optimal mean
			res = computeMuSigma(X, w, I, J);
			mu = res[0];
			sigma = res[1];
			
			int[] T = new int[k+1];
			int maxRepetitions = L-1;

			for (int i=0; i<k; i++)
			    if (J[i] < I[i]) maxRepetitions = 1;

			for (int m=1; m<=maxRepetitions; m++) {
			    T[k] = L-m-1;
			    			    
			    //initialize U indices
			    for (int i=0; i<k; i++) {
				if (J[i] == 0) {
				    U[i] = 0;
				    T[i] = 0;
				} else { 
				    U[i] = J[i]-1;
				    T[i] = J[i]-1;
				}
			    }
		    
			    double c = Double.POSITIVE_INFINITY;
		    
			    //iterate over all indices u_i in {j_i-1,j_i}
			    while (true) {
			    
				double t = C.get(T);

				if (t < c) {
				    c = t;
				    meanList = new ArrayList<ArrayList<Double>>();
				    for (ArrayList<Double> x : Z.get(T)) {
					meanList.add(new ArrayList<Double>(x));
				    }
				} else if (t == c) {
				    for (ArrayList<Double> x : Z.get(T)) {
					if (!meanList.contains(x)) {
					    meanList.add(new ArrayList<Double>(x));
					}
				    }
				}
			
				// increment U indices
				carry = true;
				d = 0;
				while (carry && (d < k)) {
				    if (U[d] == J[d]) {
					if (U[d] > 0) {
					    U[d]--;
					    T[d]--;
					}
					d++;
				    } else {
					U[d]++; 
					T[d]++;
					carry = false;
				    }
				}
				if (carry) break; //finished U iteration
			    }

			    double cost = m*sigma + c;
			    if (cost < C.get(H)) {
				C.put(H, cost);
				ArrayList<ArrayList<Double>> Q = new ArrayList<ArrayList<Double>>();
				for (ArrayList<Double> x : meanList) {
				    for (int i=0; i<m; i++) x.add(mu);
				    if (!Q.contains(x)) Q.add(x);
				}
				Z.put(H, Q);
			    } else if (cost == C.get(H)) {
				for (ArrayList<Double> x : meanList) {
				    for (int i=0; i<m; i++) x.add(mu);
				    if (!Z.get(H).contains(x))
					Z.get(H).add(new ArrayList<Double>(x));
				}
			    }
			}
		
			//increment J indices
			carry = true;
			d = 0;
			while (carry && d < k) {
			    if (J[d] == I[d]) {
				J[d] = 0;
				d++;
			    } else {
				J[d]++; 
				carry = false;
			    }
			}
			if (carry) break; //finished J iteration
		    }
		}
	    }
	}

	//return results
	Double[] Fvalues = new Double[maxLength];
	ArrayList<ArrayList<ArrayList<Double>>> means = new ArrayList<ArrayList<ArrayList<Double>>>();

	int[] r = new int[k+1];
	for (int i=0; i<k; i++) r[i] = N[i]-1;

	for (int l=0; l<maxLength; l++) {
	    r[k] = l;
	    Fvalues[l] = C.get(r);
	    means.add(Z.get(r));
	}
	return new Pair<Double[], ArrayList<ArrayList<ArrayList<Double>>>>(Fvalues, means);
    }

    
    /**
     * A helper function to compute the values mu and sigma for given indices I and J.
     */
    public static double[] computeMuSigma(double[][] X, double[] w, int[] I, int[] J) {
	double count = 0.0;
	double mu = 0.0;
	for (int i=0; i<X.length; i++) {
	    count += w[i] * (I[i] - J[i] + 1.0);
	    double temp = 0.0;
	    for (int j=J[i]; j<=I[i]; j++)
		temp += X[i][j];
	    mu += w[i] * temp;
	}
	mu = mu / count;

	double sigma = 0.0;
	for (int i=0; i<X.length; i++) {
	    double temp = 0.0;
	    for (int j=J[i]; j<=I[i]; j++)
		temp += Math.pow(X[i][j]-mu, 2);
	    sigma += w[i] * temp;
	}

	return new double[] {mu, sigma};
    }
}
