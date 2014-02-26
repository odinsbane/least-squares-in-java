package org.orangepalantir.leastsquares.fitters;

import org.orangepalantir.leastsquares.Function;
import org.orangepalantir.leastsquares.Fitter;
import Jama.LUDecomposition;
import Jama.Matrix;

public class NonLinearSolver implements Fitter {
    /**
     * Ignores second derivatives, not very good.
     * for solving non-linear least squares fit of f(x:A) = z with sets of x,z data points.
     * */
    double[][] X;    //values
    double[]   A,    //Parameter Set
               Z;    //output vaues
    
    Function FUNCTION;
    double[] ERROR;
    double[][] DERIVATIVES;
    
    double[][] HESSIAN;
    
    double DELTA = 0.000001;   //used for calculating derivatives 
    double MINERROR = 1e-9; //for evaluating 
    double MINCHANGE = 1e-3;    //minimumchanges
    public NonLinearSolver(Function funct){
        
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
    
    public void calculateDerivatives(){
        double[] abefore = new double[A.length];
        double[] aafter = new double[A.length];
        System.arraycopy(A,0,abefore,0,A.length);
        System.arraycopy(A,0,aafter,0,A.length);
        for(int j = 0; j<A.length; j++){
            aafter[j] += DELTA;
            abefore[j] -= DELTA;
            
            if(j>0){
                aafter[j-1] = A[j-1];
                abefore[j-1] = A[j-1];
            }
                
            for(int i = 0; i<Z.length; i++){
                DERIVATIVES[i][j] = (FUNCTION.evaluate(X[i],aafter) - FUNCTION.evaluate(X[i],abefore))/(2*DELTA);
            }
        }
    }
    
    
    /**
     *  Takes the current error, and the current parameter set and calculates the
     *  changes, then returns the maximum changed value
     * */
    public double iterateValues(){
        double[][] matrix = new double[A.length][A.length];
        double[] right = new double[A.length];
        for(int i = 0; i<A.length; i++){
            for(int k = 0; k<Z.length;k++){
                    right[i] += ERROR[k]*DERIVATIVES[k][i];
                    for(int l=0; l<A.length; l++){
                        matrix[i][l] += DERIVATIVES[k][i]*DERIVATIVES[k][l];
                    }
            }
                
        }
        
        
        
        Matrix coeff = new Matrix(matrix);
        Matrix sols = new Matrix(right,A.length);
        
        for(int i = 0; i<matrix.length; i++){
            for(int j=0; j<matrix[0].length; j++){
                System.out.print(matrix[i][j] + "\t");
            }
            System.out.println("|\t" + right[i]);
        }
        
        LUDecomposition adecom = coeff.lu();

        Matrix out = adecom.solve(sols);
        
        double[][] values = out.getArray();
        
        double max_change = Math.abs(values[0][0]);
        for(int i = 0; i<A.length; i++){
            max_change = max_change>Math.abs(values[i][0])?max_change:Math.abs(values[i][0]);
            A[i] += values[i][0]*0.01;
        }

        return max_change;
        
        
    }
    public void iterateValuesB(){
        Matrix system = new Matrix(DERIVATIVES);
        Matrix params = new Matrix(ERROR,ERROR.length);
        //params = params.times(-1.);
        Matrix out = system.solve(params);
        
        double[][] values = out.getArray();

        for(int i = 0; i<A.length; i++){
            A[i] += values[i][0]*0.001;
        }

    }
    
    
    public void initializeWorkspace(){
        ERROR = new double[Z.length];
        DERIVATIVES = new double[Z.length][A.length];
    }
    
    /**
     *  Main routine, call this and the Parameters are iterated until it is finished.
     * */
    public void fitData(){
        initializeWorkspace();
        int i;
        double changes = 0;
        for(i = 0; i<10000; i++){

            double e = calculateErrors();
            if(e<MINERROR)
                break;
            calculateDerivatives();
            try{
                changes = iterateValues();
                if(changes<MINCHANGE)
                    break;
            } catch(Exception exc){
                break;
            }
            System.out.println("changes: " + changes);
        }
        System.out.println("After " + i + " iterations");
        
        
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

}
