package org.orangepalantir.leastsquares.fitters;

import org.orangepalantir.leastsquares.Fitter;

import Jama.Matrix;
import org.orangepalantir.leastsquares.Function;

public class MarquardtFitter implements Fitter {
    /**
     * for solving non-linear least squares fit of f(x:A) = z with sets of x,z data points.
     * */
    double[][] X;    //values
    double[]   A,    //Parameter Set
               Z;    //output vaues

    Function FUNCTION;
    double[] ERROR;

    double[][] DERIVATIVES;
    double[][] ALPHA_PRIME;
    double[] LAMBDA;
    double[] BETA;

    double DELTA = 0.000001;   //used for calculating derivatives
    double MINERROR = 1e-9; //for evaluating
    double MINCHANGE = 1e-3;    //minimumchanges

    public int ITERATIONS = 0;
    public MarquardtFitter(Function funct){

        FUNCTION=funct;

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
     * @param k - index of the parameter that the derivative is being taken of
     * @param params - array that will be used to calculate derivatives, note params will be modified.
     * @param x - set of values to use.
     * @return derivative of function
     **/
    public double calculateDerivative(int k, double[] params, double[] x){
        double b, a;

        params[k] -= DELTA;
        b = FUNCTION.evaluate(x, params);
        params[k] += 2*DELTA;
        a = FUNCTION.evaluate(x, params);
        params[k] -= DELTA;
        return (a - b)/(2*DELTA);

    }

    /**
     *  Creates an array of derivatives since each one is used 3x's
     **/
    public void calculateDerivatives(){
        double[] working = new double[A.length];
        System.arraycopy(A,0,working,0,A.length);
        for(int j = 0; j<A.length; j++){

            for(int i = 0; i<Z.length; i++){
                DERIVATIVES[i][j] = calculateDerivative(j, working, X[i]);
            }
        }


    }

    /**
     *
     *  NOT USED.
     *  Given a set of parameters, and inputs, calculates the
     *  second derivatives
     *  d^2/d a_l d a_k (F)
     *
     * @param k - index of the parameter that the derivative is being taken of
     * @param params - array that will be used to calculate derivatives, note params will be modified.
     * @param x - set of values to use.
     **/
    public double calculateSecondDerivative(int l, int k, double[] params, double[] x){
        double b, a;

        params[l] -= DELTA;
        b = calculateDerivative(k, params, x);

        params[l] += 2*DELTA;
        a = calculateDerivative(k, params, x);


        params[l] -= DELTA;
        return (a - b)/(2*DELTA);

    }

    public void createBetaMatrix(){
        BETA = new double[A.length];
        double[] working = new double[A.length];
        System.arraycopy(A,0,working,0,A.length);

        for(int k = 0; k<BETA.length; k++){
            for(int i = 0; i<X.length; i++){

                BETA[k] += ERROR[i]*DERIVATIVES[i][k];

            }
        }

    }

    public void createAlphaPrimeMatrix(){
        ALPHA_PRIME = new double[A.length][A.length];

        int n = A.length;
        for(int k = 0; k<n; k++){
            for(int l = 0; l<n; l++){

                    for(int i = 0; i<X.length; i++)
                        ALPHA_PRIME[l][k] += DERIVATIVES[i][k] * DERIVATIVES[i][l];

                    if(k==l)
                        ALPHA_PRIME[l][k] = ALPHA_PRIME[l][k]*(1 + LAMBDA[k]);

            }

        }



    }

    /**
     *  Takes the current error, and the current parameter set and calculates the
     *  changes, then returns the maximum changed value
     * */
    public void iterateValues(){

        calculateDerivatives();

        createBetaMatrix();
        createAlphaPrimeMatrix();

        Matrix alpha_matrix = new Matrix(ALPHA_PRIME);
        Matrix beta = new Matrix(BETA, BETA.length);

        Matrix out = alpha_matrix.solve(beta);

        double[][] delta_a = out.getArray();

        for(int i = 0; i<A.length; i++)
            A[i] += delta_a[i][0];




    }



    public void initializeWorkspace(){
        ERROR = new double[Z.length];
        DERIVATIVES = new double[Z.length][A.length];
        LAMBDA = new double[A.length];
        for(int i = 0; i<A.length; i++)
            LAMBDA[i] = 0.01;
    }

    /**
     *  Main routine, call this and the Parameters are iterated until it is finished.
     * */
    public void fitData(){
        initializeWorkspace();
        double error = calculateErrors();
        double nerror, value;
        double[] acopy = new double[A.length];

        int i;
        for(i = 0; i<10000; i++){

            try{
                System.arraycopy(A,0,acopy,0,A.length);
                iterateValues();

            } catch(Exception exc){
                System.out.println("Broke after " + i + " iterations");
                printMatrix();
                //printMatrix();
                throw new RuntimeException(exc);
            }

            nerror = calculateErrors();

            value = error - nerror;


            if(value<0){

                //reject changes
                System.arraycopy(acopy, 0, A, 0, acopy.length);
                updateLambda(value);
                nerror = calculateErrors();
            } else{


                if(value<0.0001){
                    ITERATIONS = i;
                    break;
                }
            }

            error = nerror;
        }


    }

    public void updateLambda(double value){

        if(value<0){
            for(int i = 0; i<LAMBDA.length; i++)
                LAMBDA[i] = LAMBDA[i]*10;


        } else {

            for(int i = 0; i<LAMBDA.length; i++)
                LAMBDA[i] = LAMBDA[i]*0.1;

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
        return new double[0];
    }

    /**
     * for stack traces
     */
    public String printMatrix(){
        String message = "";
        for(int i = 0; i<ALPHA_PRIME.length; i++){
            for(int j = 0; j<ALPHA_PRIME[0].length; j++){
                message += ALPHA_PRIME[i][j] + "\t";

            }

            message += "| " + BETA[i] + "\n";
        }
        return message;
    }

}