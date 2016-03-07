package org.orangepalantir.leastsquares.examples;

import org.orangepalantir.leastsquares.Fitter;
import org.orangepalantir.leastsquares.Function;
import org.orangepalantir.leastsquares.fitters.LinearFitter;

import java.util.Arrays;

/**
 * Created by msmith on 3/7/16.
 */
public class ParabolaExample {

    public static void main(String[] args){
        double[][] xs = {
                {0}, {1}, {2}, {3}, {4}, {5}
        };
        double[] zs = {1.0, 0.9, 1.0, 1.3, 1.8, 2.5};

        Function fun = new Function(){
            @Override
            public double evaluate(double[] values, double[] parameters) {
                double A = parameters[0];
                double B = parameters[1];
                double C = parameters[2];
                double x = values[0];
                return A*x*x + B*x + C;
            }
            @Override
            public int getNParameters() {
                return 3;
            }

            @Override
            public int getNInputs() {
                return 1;
            }
        };

        Fitter fit = new LinearFitter(fun);
        fit.setData(xs, zs);
        fit.setParameters(new double[]{0,0,0});

        fit.fitData();

        System.out.println(Arrays.toString(fit.getParameters()));
    }

}
