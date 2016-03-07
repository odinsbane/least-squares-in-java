package org.orangepalantir.leastsquares.fitters;

import org.junit.Test;
import org.junit.Assert;
import org.orangepalantir.leastsquares.functions.LinearFunction;

/**
 *
 *
 * Created by msmith on 2/26/14.
 */
public class LinearFitterTest {
    @Test
    public void testLinearFitter(){
        int N = 100;
        double[] z = new double[N];
        double[][] x = new double[N][2];
        for(int i = 0; i<N; i++){
            double v = 0.1*(N/2 - i);
            z[i] = 1 - 2* v + v*v;
            x[i][0] = v;
            x[i][1] = v*v;
        }
        LinearFitter fit = new LinearFitter(new LinearFunction(2));
        fit.setData(x, z);
        fit.setParameters(new double[]{1,1,1});
        fit.fitData();
        double[] output = fit.getParameters();

        Assert.assertEquals(1.0, output[0], 1.0e-15);
        Assert.assertEquals(-2.0, output[1], 1.0e-14);
        Assert.assertEquals(1.0, output[2], 1.0e-14);

    }

}