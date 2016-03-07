least-squares-in-java
=====================

Java Least Squares fitting library.
-----------------------------------

This is a small least squares fitting library made in java.
It was originally used in the development of an image analysis tool 
SpeckleTrackerJ. http://athena.physics.lehigh.edu/speckletrackerj/ It uses
methods described in "Numerical Recipes in C", Second Edition (1992).

The outline of using this library is you have a set of data that you
would like fit a function to. Define the Function, how it evaluates
the input and paramters. Create a Fitter (three types currently) initialized
to the function. Set the fitter with the data and make an initial guess of the
parameters. Tell the fitter to find the best set of parameters which minimizes the
sum of squares of error.

Example
-------

Lets say we have some data.

    double[][] xs = {
        {0}, {1}, {2}, {3}, {4}, {5}
    };
    double[] z = {1.1, 0.9, 1.0, 1.35, 1.82, 2.5};

And we want to fit a quadratic function f(x) = A*x*x + B*x + C. The parameters are A, B, C and the
input is x. So we have 3 parameters and 1 input.

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

The next step is to create a fitter.

    Fitter fit = new LinearFitter(fun);

Set the data, and make a guess at the initial parameters.

    fit.setData(xs, zs);
    fit.setParameters(new double[]{0,0,0});

    fit.fitData();

The results can be obtained by using getParameters.

    System.out.println(Arrays.toString(fit.getParameters()));
    #[0.10000000000000023, -0.20000000000000123, 1.000000000000001]

This example is included in the examples package.


Fitter implementations
----------------------

LinearFitter: for use with equations that are linearly dependent on the input parameters. Standard
linear regression where derivatives are taken by setting all of the parameters to zero, except for the
 one of interest.

 NonLinearSolver: Similar to the LinearSolver, except it will run multiple iterations, there is a damping
  factor, and the derivatives are calculated by taken by varying the parameters a small amount from the
  previous guess.

 MarquardtFitter: Similar to the NonLinearSolver except the damping is adaptive.
