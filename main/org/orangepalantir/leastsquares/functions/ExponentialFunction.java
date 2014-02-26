package org.orangepalantir.leastsquares.functions;

import org.orangepalantir.leastsquares.Function;

/**
 *
 * Created by msmith on 2/26/14.
 */
public class ExponentialFunction implements Function {
    @Override
    public double evaluate(double[] values, double[] parameters) {
        return parameters[0] + parameters[1]*Math.exp(parameters[2]*values[0]);
    }

    @Override
    public int getNParameters() {
        return 3;
    }

    @Override
    public int getNInputs() {
        return 1;
    }
}
