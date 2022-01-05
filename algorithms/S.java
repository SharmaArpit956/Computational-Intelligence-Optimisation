package algorithms;

import static utils.algorithms.Misc.*;
import interfaces.Algorithm;
import interfaces.Problem;
import utils.RunAndStore.FTrend;

public class S extends Algorithm

{
	@Override
	public FTrend execute(Problem problem, int maxEvaluations) throws Exception {

		// we always need an object of the
		// kind FTrend (for storing the fitness trend), and variables
		// for storing the dimesionality value and the bounds
		// of the problem as shown below
		FTrend FT = new FTrend();

		// variable for storing the number of dimensions of the problem
		int problemDimension = problem.getDimension();

		// variable for storing the boundaries of each component of the solution
		double[][] bounds = problem.getBounds();

		double[] x = new double[problemDimension];// variable for storing the solution

		// variable for storing the best solution available
		double[] xBest = new double[problemDimension];

		// Fixed 5 variable solution for checking the deterministic nature of the
		// algorithm (fitness value should be same on multiple runs)
		// double[] xBest = { 1, 1, 1, 1, 1 };

		double fBest; // variable for storing the best fitness value available

		int i = 0; // The counter for keeeping track of the computational budget consumed

		// if the initial solution, we are getting from somwehere else like, from
		// another optimiser
		if (initialSolution != null) {

			// the best solution till now is the first solution as we donot have any other
			// solution yet
			xBest = initialSolution;

			// Similarly, the best fitness value till now is the first finess value as we
			// donot any other finess value yet
			fBest = initialFitness;

		} else// random intitial guess
		{
			// generate Random Solution within the search
			// space and as it is the only solution till
			// now, it will be the best solution we know
			xBest = generateRandomSolution(bounds, problemDimension);

			// calculate the fitness value of the best solution available till now.
			fBest = problem.f(xBest);

			FT.add(i, fBest);// store the initital guess

			i++; // increase the counter after every fitness functional call
		}

		// temp variables
		double[] radius = new double[problemDimension]; // exploratory radius
		double alpha = 0.4; // for calculating exploratory radius

		// makes exploratory radius equal to α(Upper limit -
		// lower limit) for every component of the solution
		for (int p = 0; p < problemDimension; p++) {
			radius[p] = alpha * (bounds[p][1] - bounds[p][0]);
		}

		double fx = fBest;// initialize fitness function
		x = xBest; // initial solution

		// Flag for determining if the fitness ever improved or not. As no fitness value
		// improvement has been made till now, it is initialised as true
		boolean neverImproved = true;

		// S Algorithm here
		while (i < maxEvaluations) { // keep the loop running till the budget expires

			neverImproved = true; // re-initialize the flag to true for the next iteration

			// for every axis of the variable and till the budget expires
			for (int j = 0; j < problemDimension && i < maxEvaluations; j++) {

				// perturb the jth component (starting from 0) towards negative jth axis
				// treating other components as constants (like partial derivative)
				x[j] = xBest[j] - radius[j];

				// After every perturbation of the solution,a correction method is needed as as
				// there are chances of the component of having landed outside the search space
				x = toro(x, bounds);

				// fitness value calculaion acoording to the problem being handled
				fx = problem.f(x);
				i++; // Increase the counter after every fitness functional call

				// Fitness Improved (Found a lower point on the curve for the jth dimension)
				if (fx <= fBest) {

					// Very important. this means that improvement in fitness value happened.
					// Therefore, we cannot say that improvement never happened
					neverImproved = false;

					xBest[j] = x[j]; // As the improvenment happened, the best solution is updated
					fBest = fx; // Similarly, the best fitness value is also updated

					// the new best fitness value is added to the fitness trend array
					FT.add(i, fBest);

				} else { // when the fitness value doesnt improve

					// it symbolises coming back to the same position if fitness value
					// doesnt improve in the pseudo code, so that it doesnt get worse than what it
					// initially was atleast.
					// but, it means nothing in programmming
					// x[j] = xBest[j];

					// As the fitness value didnt improve, move half the exploratory
					// radius in the opposite direcion
					x[j] = xBest[j] + (radius[j] / 2);

					// After every perturbation of the solution, correction method is needed as as
					// there are chances of x landing outside the search space in both cases
					x = toro(x, bounds);

					fx = problem.f(x); // fitness functional call
					i++; // Increase the counter after every fitness functional call

					// Fitness Improved (Found a lower point on the curve for the jth dimension)
					if (fx <= fBest) {

						// Very important. this means that improvement in fitness value happened.
						// Therefore, we cannot say improvement never happened
						neverImproved = false;

						xBest[j] = x[j]; // As the improvenment happened, the best solution is updated
						fBest = fx; // Similarly, the best fitness value is also updated

						// the new best fitness value is added to the fitness trend array
						FT.add(i, fBest);

					} else { // if the fitness value doesnt improve
						x[j] = xBest[j]; // store the best solution as the current solution

						// store the best fitness functions as the current fitness function
						fx = fBest;
					}
				}
			}

			// If the fitness function never improved throughout the for loop
			if (neverImproved == true) {

				for (int d = 0; d < problemDimension; d++) {

					// half the radius to do the perturbations as it means there
					// are more chances of the minima lying near the current point rather
					// than far from it(which we have checked and not found at this point)
					radius[d] = radius[d] / 2;
				}
			}
		}

		finalBest = xBest; // save the final best
		FT.add(i, fBest); // the best fitness value is added to the fitness trend array

		return FT; // return the fitness trend
	}

}
