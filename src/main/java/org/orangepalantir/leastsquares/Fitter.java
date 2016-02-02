package org.orangepalantir.leastsquares;

public interface Fitter {
    
    public void setData(double[][] xvalues, double[] zvalues);
    public void setParameters(double[] parameters);
    public double[] getParameters();
    public double[] getUncertainty();
    public void fitData();
    public double calculateErrors();
}
