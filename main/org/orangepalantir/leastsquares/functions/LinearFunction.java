package org.orangepalantir.leastsquares.functions;

import org.orangepalantir.leastsquares.Function;

/**
 * Generic linear function, with an offset.
 *
 * A + B x1 + C x2 + D x3 ...
 *
 * Where the number of inputs, x's is indicated in the constructor and
 * the number of parameters is the number of terms plus 1.
 *
 * Created by msmith on 2/26/14.
 */
public class LinearFunction implements Function {
    int terms;
    public LinearFunction(int terms){
        this.terms = terms;
    }
    @Override
    public double evaluate(double[] values, double[] parameters) {
        double sum = parameters[0];
        for(int i = 0; i<values.length; i++){
            sum += parameters[i+1]*values[i];
        }
        return sum;
    }

    @Override
    public int getNParameters() {
        return terms + 1;
    }

    @Override
    public int getNInputs() {
        return terms;
    }
}
