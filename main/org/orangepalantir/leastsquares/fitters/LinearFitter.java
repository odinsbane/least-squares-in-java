package org.orangepalantir.leastsquares.fitters;

import org.orangepalantir.leastsquares.Fitter;
import Jama.Matrix;
import org.orangepalantir.leastsquares.Function;

/**
 * for solving linear least squares fit of f(x:A) = z with sets of x,z data points.
  * derivatives are evaluated numerically and the function is assumed to be linear in A.
  * Then a linear least squares fit is performed.
 * */
public class LinearFitter implements Fitter {
    private final double[] working;

    double[][] X;    //values
    double[]   A,    //Parameter Set
               Z;    //output vaues

    Function FUNCTION;
    double[] ERROR;

    double[][] DERIVATIVES;
    double[] BETA;
    double[][] ALPHA;
    double DELTA = 0.000001;   //used for calculating derivatives

    public LinearFitter(Function funct){
        FUNCTION=funct;
        working = new double[funct.getNParameters()];
    }


    /**
     *      Sets the values of of the original data points that are going to be fit.
     * */
    public void setData(double[][] xvalues, double[] zvalues){

        if(xvalues.length != zvalues.length)
            throw new IllegalArgumentException("there must be 1 z value for each set of x values");
        else if(xvalues[0].length != FUNCTION.getNInputs())
            throw new IllegalArgumentException("The length of parameters is longer that the parameters accepted by the function");
        X = xvalues;
        Z = zvalues;

    }

    /**
     *
     * */
    public void setParameters(double[] parameters){
        if(parameters.length != FUNCTION.getNParameters())
            throw new IllegalArgumentException("the number of parameters must equal the required number for the function: " + FUNCTION.getNParameters());
        A = new double[parameters.length];
        System.arraycopy(parameters,0,A,0,parameters.length);
    }

     /**
      *     returns the cumulative error for current values
      **/
    public double calculateErrors(){
        double new_error = 0;
        for(int i = 0; i<Z.length; i++){

            double v = FUNCTION.evaluate(X[i],A);
            ERROR[i] = Z[i] - v;
            new_error += Math.pow(ERROR[i],2);
        }
        return new_error;

    }

    /**
     *  Given a set of parameters, and inputs, calculates the derivative
     *  of the k'th parameter
     *  d/d a_k (F)
     *
     *  The function is linear in parameters eg:
     *  f(x) = params[0]* f0(x) + params[1]*f1(x) + ...
     *
     *  So the derivative k is fk(x) or the function evaluated when all of the parameters
     *  are zero except for the k'th parameter which is 1.
     *
     *
     * @param k - index of the parameter that the derivative is being taken of
     * @param x - set of values to use.
     * @return derivative of function
     **/
    public double calculateDerivative(int k, double[] x){
        for(int i = 0; i<FUNCTION.getNParameters(); i++){
            working[i] = i==k?1:0;
        }
        return FUNCTION.evaluate(x, working);

    }

    /**
     *  Creates an array of derivatives since each one is used 3x's
     **/
    public void calculateDerivatives(){
        for(int j = 0; j<A.length; j++){
            for(int i = 0; i<Z.length; i++){
                DERIVATIVES[i][j] = calculateDerivative(j, X[i]);
            }
        }
    }


    public void createBetaMatrix(){
        BETA = new double[A.length];

        for(int k = 0; k<BETA.length; k++){
            for(int i = 0; i<X.length; i++){

                BETA[k] += ERROR[i]*DERIVATIVES[i][k];

            }
        }

    }

    public void createAlphaMatrix(){
        ALPHA = new double[A.length][A.length];

        int n = A.length;
        for(int k = 0; k<n; k++){
            for(int l = 0; l<n; l++){

                    for(int i = 0; i<X.length; i++)
                        ALPHA[l][k] += DERIVATIVES[i][k] * DERIVATIVES[i][l];

            }

        }



    }

    /**
     *  Takes the current error, and the current parameter set and calculates the
     *  changes, then returns the maximum changed value
     * */
    public void iterateValues(){
        calculateErrors();
        calculateDerivatives();

        createBetaMatrix();
        createAlphaMatrix();

        Matrix alpha_matrix = new Matrix(ALPHA);
        Matrix beta = new Matrix(BETA, BETA.length);
        
        Matrix out = alpha_matrix.solve(beta);

        double[][] delta_a = out.getArray();

        for(int i = 0; i<A.length; i++)
            A[i] += delta_a[i][0];

    }

    public void printMatrix(){
        for(int i = 0; i<ALPHA.length; i++){
            for(int j = 0; j<ALPHA[0].length; j++){
                System.out.print(ALPHA[i][j] + "\t");

            }

            System.out.println("| " + BETA[i] );
        }
    }

    public void initializeWorkspace(){
        ERROR = new double[Z.length];
        DERIVATIVES = new double[Z.length][A.length];
    }

    /**
     *  Main routine, call this and the Parameters are iterated one time because
     *  the equations are assumed to be linear in parameter space.
     * */
    public void fitData(){
        initializeWorkspace();

        try{
          iterateValues();

        } catch(Exception exc){
            printMatrix();
            exc.printStackTrace();

         }

    }


    /**
     *      Gets the current set of parameters values.
     * */
    public double[] getParameters(){
        return A;
    }


    @Override
    public double[] getUncertainty() {
        printMatrix();
        double sum_a = 0;
        double sum_b = 0;
        for(int i = 0;i<DERIVATIVES.length;i++){
            sum_a += Math.pow(DERIVATIVES[i][0],2);
            sum_b += Math.pow(DERIVATIVES[i][1],2);
        }
        Matrix a_matrix = new Matrix(ALPHA);
        Matrix b = a_matrix.inverse();

        for(int i = 0; i<b.getColumnDimension(); i++){
            for(int j = 0; j<b.getRowDimension(); j++){
                System.out.print(b.get(i,j) + "\t");
            }
            System.out.println();
        }
        System.out.println(sum_a/sum_b);
        double[] residuals = new double[A.length];
        double error = calculateErrors();
        for(int i = 0; i<A.length; i++){
            residuals[i] = error*Math.sqrt(b.get(i,i))/2;
        }
        return residuals;
    }

}