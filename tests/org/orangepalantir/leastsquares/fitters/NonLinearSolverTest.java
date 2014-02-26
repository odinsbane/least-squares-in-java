package org.orangepalantir.leastsquares.fitters;

import org.junit.Assert;
import org.junit.Test;
import org.orangepalantir.leastsquares.Fitter;
import org.orangepalantir.leastsquares.Function;
import org.orangepalantir.leastsquares.functions.ExponentialFunction;

import java.util.Arrays;

/**
 * Created by msmith on 2/26/14.
 */
public class NonLinearSolverTest {
    @Test
    public void testExponentialFunction(){
        int N = 100;
        double[] z = new double[N];
        double[][] x = new double[N][1];

        for(int i = 0; i<N; i++){
            x[i][0] = 0.1*i;
            z[i] = Math.exp(-x[i][0] / 2.0);
        }

        Function exponential = new ExponentialFunction();
        Fitter non_linear = new NonLinearSolver(exponential);
        non_linear.setData(x,z);
        non_linear.setParameters(new double[]{0,1,-1});
        non_linear.fitData();

        double[] results = non_linear.getParameters();
        System.out.println(Arrays.toString(results));
        Assert.assertEquals(0, results[0], 1e-12);
        Assert.assertEquals(1, results[1], 1e-12);
        Assert.assertEquals(-0.5, results[2], 1e-12);
    }
}
