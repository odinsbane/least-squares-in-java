package org.orangepalantir.filters;
import Jama.LUDecomposition;
import Jama.Matrix;

/** Class implements the Savitzy-Golay filter.  The algorithm is based on
 * the algorithm used in 
 *
 * Numerical Methods in C: The art of scientific computing second edition
 *
 * This filter is used for smoothing data.  The current implementation is
 * a 1-d approach and a symmetric 2-d approach.
 *
 *
 *     @Author: Matthew B. Smith
 *     @Url: http://orangepalantir.org
 *     @Version: 0.8
 *     @Date: 2/10/2010
 *
 **/
public class SavitzyGolayFilter{
    double[] coefficients;
    int LEFT, RIGHT, ORDER;
    /**
     * Prepares the filter coefficients
     *
     * @param left number of places to the left of the origion
     * @param right number of places to the right of the origion
     * @param order of polynomial used to fit data
     * */
    public SavitzyGolayFilter(int left, int right, int order){

        coefficients = SavitzkyGolayCoefficients(left,right,order);
        LEFT = left;
        RIGHT = right;
        ORDER = order;
    }

    /**
     * Prepares the filter coefficients
     *
     * @param left number of places to the left of the origion
     * @param right number of places to the right of the origion
     * @param order of polynomial used to fit data
     * @param derivative use a SavitzkyGolayFilter to calculate the derivative or order.
     * */
    public SavitzyGolayFilter(int left, int right, int order, int derivative){

        coefficients = SavitzkyGolayCoefficients(left,right,order, derivative);
        LEFT = left;
        RIGHT = right;
        ORDER = order;
    }


    /**
     * filters the data by assuming the ends are reflected.
     *
     * @param data the array that is going to be filtered.
     * */
    public double[] filterData(double[] data){
        int pts = data.length;
        double[] ret_value = new double[pts];

        int j;
        //reflected
        for(j = 0; j<LEFT; j++)
            ret_value[j] = rFilter(data, coefficients, j, LEFT);


        //normal
        for(j = LEFT; j<pts - RIGHT; j++)
            ret_value[j] = filter(data, coefficients, j, LEFT);


        //reflected
        for(j = pts - RIGHT; j<pts; j++)
            ret_value[j] = rFilter(data, coefficients, j, LEFT);

        return ret_value;
    }
    /**
     * filters the data by assuming the ends are reflected.
     *
     * @param data
     * */
    public float[] filterData(float[] data){
        int pts = data.length;
        float[] ret_value = new float[pts];

        int j;
        //reflected
        for(j = 0; j<LEFT; j++)
            ret_value[j] = (float)rFilter(data, coefficients, j, LEFT);


        //normal
        for(j = LEFT; j<pts - RIGHT; j++)
            ret_value[j] = filter(data, coefficients, j, LEFT);


        //reflected
        for(j = pts - RIGHT; j<pts; j++)
            ret_value[j] = (float)rFilter(data, coefficients, j, LEFT);

        return ret_value;
    }

    /**
     *  Calculates one point in the data, similar to one step of a
     *  convolution.
     *
     * @param o - data that will be convolved
     * @param mask - kernel
     * @param start - position of the 'middle point' on the data
     * @param middle - postion of the center of the kernel on the kernel
     *
     * */
    public float filter(float[] o, double[] mask, int start, int middle){
        float out = 0;
        for(int i = 0; i<mask.length; i++)
            out += (float)mask[i]*o[start - middle + i];

        return out;

    }

    /**
     *  Calculates one point in the data, similar to one step of a
     *  convolution.
     *
     * @param o - data that will be convolved
     * @param mask - kernel
     * @param start - position of the 'middle point' on the data
     * @param middle - postion of the center of the kernel on the kernel
     *
     * */
    public double filter(double[] o, double[] mask, int start, int middle){
        double out = 0;
        for(int i = 0; i<mask.length; i++)
            out += mask[i]*o[start - middle + i];

        return out;

    }

    /**
     *  Calculates one point in the data, similar to one step of a
     *  convolution.  This one will wrap the ends of the mask.  Something
     *  like a convolution with reflected boundary conditions.
     *
     * @param o - data that will be convolved
     * @param mask - kernel
     * @param start - position of the 'middle point' on the data
     * @param middle - postion of the center of the kernel on the kernel
     *
     * */
    public double rFilter(double[] o, double[] mask, int start, int middle){
        double out = 0;
        int dex;
        for(int i = 0; i<mask.length; i++){
            dex = Math.abs(start - middle + 1 +i);
            dex = (dex<o.length)?dex: 2*o.length - dex - 1;
            out += mask[i]*o[dex];
        }

        return out;

    }

    /**
     *  Calculates one point in the data, similar to one step of a
     *  convolution.  This one will wrap the ends of the mask.  Something
     *  like a convolution with reflected boundary conditions.
     *
     * @param o - data that will be convolved
     * @param mask - kernel
     * @param start - position of the 'middle point' on the data
     * @param middle - postion of the center of the kernel on the kernel
     *
     * */
    public float rFilter(float[] o, double[] mask, int start, int middle){
        double out = 0;
        int dex;
        for(int i = 0; i<mask.length; i++){
            dex = Math.abs(start - middle + 1 +i);
            dex = (dex<o.length)?dex: 2*o.length - dex - 1;
            out += mask[i]*o[dex];
        }

        return (float)out;

    }

    /**
     *
     *   Creates the array of coefficients for using a least squares polynomial fit
     *   of the data.  Due to the fact the points are equally spaced
     *   and the polynomial least squares is linear in the coefficients,
     *   a set of coefficients can be used for every element.
     *
     *   @param nl - spaces to the left
     *   @param nr - spaces to the right
     *   @param order - the order of the polynomial being used, not that nl+nr+1 must be more than the order of the polynomial
     *   @return a 'time kernel' for convolving in time.
     * */


    public static double[] SavitzkyGolayCoefficients(int nl, int nr, int order){
        if( nl + nr + 1<= order)
            throw new java.lang.IllegalArgumentException(" The order of polynomial cannot exceed the number of points being used." +
                    "If they are equal the input equals the output.");
        int N = nr + nl + 1;
        double[] ret_values = new double[N];
        double[] xvalues = new double[N];
        for(int i = 0; i<N; i++){

            xvalues[i] = -nl + i;

        }

        int counts = 2*order+1;
        double[] moments = new double[counts];
        for(int i = 0; i<counts; i++){
            for(int j = 0; j<N; j++){

                moments[i] += Math.pow(xvalues[j],i);

            }

            moments[i] = moments[i]/N;
        }



        double[][] matrix = new double[order+1][order+1];

        for(int i = 0; i<order+1; i++){
            for(int j = 0; j<order+1; j++){
                matrix[i][j] = moments[counts - i - j - 1];
            }
            System.out.println("");
        }

        Matrix A = new Matrix(matrix);

        LUDecomposition lu = A.lu();

        Matrix x = new Matrix(new double[order+1],order+1);
        Matrix y;
        double[] polynomial;

        for(int i = 0; i<N; i++){

            for(int j = 0; j<order+1; j++)
                x.set( j , 0, Math.pow(xvalues[i],order - j));

            y = lu.solve(x);

            polynomial = y.getColumnPackedCopy();
            ret_values[i] = evaluatePolynomial(polynomial, xvalues[nl])/N;
        }


        return ret_values;

    }

    /**
     *
     *   Creates the array of coefficients for using a least squares polynomial fit
     *   of the data.  Due to the fact the points are equally spaced
     *   and the polynomial least squares is linear in the coefficients,
     *   a set of coefficients can be used for every element.
     *
     *   @param nl - spaces to the left
     *   @param nr - spaces to the right
     *   @param order - the order of the polynomial being used, not that nl+nr+1 must be more than the order of the polynomial
     *   @param derivative - returns a set of derivative coefficients
     * */


    public static double[] SavitzkyGolayCoefficients(int nl, int nr, int order, int derivative){
        if( nl + nr + 1<= order)
            throw new java.lang.IllegalArgumentException(" The order of polynomial cannot exceed the number of points being used." +
                    "If they are equal the input equals the output.");
        int N = nr + nl + 1;
        double[] ret_values = new double[N];
        double[] xvalues = new double[N];
        for(int i = 0; i<N; i++){

            xvalues[i] = -nl + i;

        }

        int counts = 2*order+1;
        double[] moments = new double[counts];
        for(int i = 0; i<counts; i++){
            for(int j = 0; j<N; j++){

                moments[i] += Math.pow(xvalues[j],i);

            }

            moments[i] = moments[i]/N;
        }



        double[][] matrix = new double[order+1][order+1];

        for(int i = 0; i<order+1; i++){
            for(int j = 0; j<order+1; j++){
                matrix[i][j] = moments[counts - i - j - 1];
            }
            System.out.println("");
        }

        Matrix A = new Matrix(matrix);

        LUDecomposition lu = A.lu();

        Matrix x = new Matrix(new double[order+1],order+1);
        Matrix y;
        double[] polynomial;

        for(int i = 0; i<N; i++){

            for(int j = 0; j<order+1; j++)
                x.set( j , 0, Math.pow(xvalues[i],order - j));

            y = lu.solve(x);

            polynomial = y.getColumnPackedCopy();
            ret_values[i] = evaluatePolynomial(polynomial, xvalues[nl], derivative)/N;
        }


        return ret_values;

    }
    /**
     * Creates a 2-d Kernel for using a savitzky-golay filter in 2-d.
     * It fits a 2-d polynomial to the data set if the order was 2.
     *
     * f(x) = a*x**2 + b*x + c*x*y + d*y + e*y**2 + f
     *
     * It is solved by least squares fiting '1' functions to polynomials
     * to 2-d polynomials, then these 1 functions can be summed as a linear
     * sum to create a 2d kernel.
     *
     * The kernel is a symmetric square kernel
     *
     *
     * @param size is the halfwidth of the kernel.
     * @param order is the order of polynomial being fit.
     * */
    public static double[][] getKernel(int size, int order){
        int N = 2*size+1;
        double[][] xyvalues = new double[N][2];

        for(int i = 0; i<N; i++){
            xyvalues[i][0] = i - size;
            xyvalues[i][1] = i - size;
        }
        int M = 0;
        for(int i = 1; i<=order+1; i++)
            M += i;

        double[][] matrix = new double[M][M];
        int pqdex = 0;
        int mndex = 0;
        //the coefficients for a_(pq) for the _XX_ equation
        for(int p = 0; p<=order; p++){
            for(int q = 0; q <= (order - p); q++){

                mndex = 0;
                for(int m = 0; m<=order; m++){
                    for(int n = 0; n<=order-m; n++){

                        double sum = 0;
                        for(int i = 0; i<N; i++){
                            for(int j = 0; j<N; j++){
                                sum += Math.pow(xyvalues[i][0], p + m)*Math.pow(xyvalues[j][1], q + n);

                            }

                        }

                        matrix[pqdex][mndex] = sum;
                        mndex++;
                    }
                }
                pqdex++;
            }
        }

        Matrix mat = new Matrix(matrix);
        LUDecomposition lud = mat.lu();


        //now we want to solve all of the different functions to get our kernel
        double[][] kernel = new double[N][N];
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){

                double[] equals = new double[M];

                int xp = 0;
                int yp = 0;

                for(int k = 0; k<M; k++){

                    equals[k] = Math.pow(xyvalues[i][0], xp)*Math.pow(xyvalues[j][1],yp);
                    yp++;
                    if(yp>order - xp){

                        yp = 0;
                        xp++;

                    }
                }



                Matrix rs = new Matrix(equals, M);
                Matrix ls = lud.solve(rs);


                //evaluate 2-d function;
                double sum = 0;

                xp = 0;
                yp = 0;

                for(int k = 0; k<M; k++){

                    sum += ls.get(k,0)*Math.pow(xyvalues[size][0], xp)*Math.pow(xyvalues[size][1],yp);


                    yp++;
                    if(yp>order - xp){

                        yp = 0;
                        xp++;

                    }
                }

                kernel[j][i] = sum;
            }

        }

        return kernel;
    }

    /**
     *      Evaluates a polynomial where:
     *      p(x) = poly[0] * x**m + poly[1] * x**(m-1) + ... + poly[m]
     *
     *  @param poly - double array representation of a polynomial
     *  @param x - the variable that will be evaluated.
     **/
    public static double evaluatePolynomial(double[] poly, double x){

        double val = 0;
        int m = poly.length;
        for(int j = 0; j<m; j++){
            val += Math.pow(x,m-j-1)*poly[j];
        }
        return val;
    }

    /**
     *      Evaluates derivative of polynomial where with:
     *      p(x) = poly[0] * x**m + poly[1] * x**(m-1) + ... + poly[m]
     *
     *  @param poly - double array representation of a polynomial
     *  @param x - the variable that will be evaluated.
     *  @param order - the order of the derivative. (n in description)
     *  @return m*(m-1)*...(m - n)*poly[0]*x**(m-n) + (m-1)*...(m - n - 1)*poly[1] * x**(m-2) ... + n!*poly[m - n]
     **/
    public static double evaluatePolynomial(double[] poly, double x, int order){

        double val = 0;
        int m = poly.length;
        for(int j = 0; j<m-order; j++){
            int pow = m - j - 1; //original power before derivatives
            val += Math.pow(x,pow - order)*poly[j]*pfact(pow,order);
        }
        return val;
    }

    /**
     * Partial factorial,
     * @param m power of x
     * @param n number of derivatives
     *
     * @return m*(m-1)*...(m - n)
     */
    public static double pfact(int m, int n){
        int p = 1;
        for(int i = 0; i<n; i++){
            p*=(m - i);
        }
        return p;
    }
}

