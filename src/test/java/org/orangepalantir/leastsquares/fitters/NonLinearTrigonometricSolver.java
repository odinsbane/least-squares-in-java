package org.orangepalantir.leastsquares.fitters;

import org.junit.Assert;
import org.junit.Test;
import org.orangepalantir.leastsquares.Function;

public class NonLinearTrigonometricSolver {

  // Solves the following non-linear set of equations:

  //  a - e + bcosθ1 + csinθ1 +  d * sin(θ1 + θ2) ) = z

  //  ( bsinθ1 + ccosθ1 + d * cos(θ1 + θ2) ) * sinθ0 = x

  //  ( bsinθ1 + ccosθ1 + d * cos(θ1 + θ2) ) * cosθ0 = y

  // given x, y, z, solve for θ0, θ1, θ2

  static final double a = 125;
  static final double b = 143;
  static final double c = 50;
  static final double d = 142;
  static final double e = 96;

  static final double x = 0;
  static final double y = 192;
  static final double z = 172;

  @Test
  public void testNonLinearTrigonometricSolver() {

    double[][] xs = { { -1 }, { 0 }, { 1 } };

    double[] zs = { z, x, y };

    double r = Math.sqrt(x * x + y * y);
    final double sinTheta0 = x / r;
    final double cosTheta0 = y / r;

    Function f = new Function() {

      @Override
      public double evaluate(double[] values, double[] parameters) {

        double t1 = parameters[0];
        double t2 = parameters[1];
        if (values[0] == -1) {
          return a - e + b * Math.cos(t1) + c * Math.sin(t1) + d * Math.sin(t2 + t1);

        } else if (values[0] == 0) {
          return (b * Math.sin(t1) + c * Math.cos(t1) + d * Math.cos(t2 + t1)) * sinTheta0;
        } else {
          return (b * Math.sin(t1) + c * Math.cos(t1) + d * Math.cos(t2 + t1)) * cosTheta0;
        }
      }

      @Override
      public int getNParameters() {

        return 2;
      }

      @Override
      public int getNInputs() {

        return 1;
      }
    };

    NonLinearSolver fit = new NonLinearSolver(f);
    fit.setData(xs, zs);
    double[] params = { 0, 0 };
    fit.setParameters(params);

    fit.fitData();
    // improving results.
    fit.setMinChange(1e-32);
    fit.setMinError(1e-32);
    fit.setStepSize(0.5);

    fit.fitData();

    double t1 = fit.getParameters()[0];
    double t2 = fit.getParameters()[1];
    double arg = y / (b * Math.sin(t1) + c * Math.cos(t1) + d * Math.cos(t2 + t1));
    // System.out.println(" " + arg);
    double theta0 = Math.acos(arg) * Math.signum(x);
    System.out.println(Math.toDegrees(theta0));
    System.out.println(Math.toDegrees(fit.getParameters()[0]));
    System.out.println(Math.toDegrees(fit.getParameters()[1]));

    Assert.assertEquals(0, Math.toDegrees(theta0), 1e-16);
    Assert.assertEquals(0, Math.toDegrees(fit.getParameters()[0]), 1e-16);
    Assert.assertEquals(0, Math.toDegrees(fit.getParameters()[1]), 1e-16);

  }
}