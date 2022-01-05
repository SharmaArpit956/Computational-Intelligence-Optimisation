package algorithms;

import static utils.algorithms.Misc.generateRandomSolution;
import interfaces.Algorithm;
import interfaces.Problem;
import utils.RunAndStore.FTrend;
import static utils.algorithms.Misc.*;
import utils.random.RandUtils;


public class DE extends Algorithm {

	@Override
	public FTrend execute(Problem problem, int maxEvaluations) throws Exception {

		// Scaling Factor for the difference vector for the mutatioon
		double F = getParameter("p0").doubleValue();

		double CR = getParameter("p1").doubleValue(); // Crossover Rate for calculating the crossover

		FTrend FT = new FTrend();
		int problemDimension = problem.getDimension(); // Number of dimensions of the problem

		// variable 'M' in the Pseudo code (number of individuals in the population)
		int populationSize = problemDimension * 10;

		// An array containing boundaries for each dimension
		double[][] bounds = problem.getBounds();

		// An 2d Array containing population of individuals having variables in every
		// dimension of the problem
		double[][] population = new double[populationSize][problemDimension];

		// Initialize the varible for the fitness values of each invidual of the
		// population of the corresponding
		// dimension
		double[] fitnesses = new double[populationSize];

		// initialize the fittest individual
		double[] best = new double[problemDimension];

		// initialize the best fitness value
		double fBest = Double.NaN;

		int i = 0;

		// PSEUDO CODE: Pop_g <-- randomly sample M n-dimentional individuals within D
		// Generate random individuals for the entire population within the search
		// space
		for (int j = 0; j < populationSize; j++) {

			double[] tmp = generateRandomSolution(bounds, problemDimension);
			for (int n = 0; n < problemDimension; n++)
				population[j][n] = tmp[n];

			// The fitness value of jth invidual of the population
			fitnesses[j] = problem.f(population[j]);

			// if its the first individual or the fitness value of the jth individual is
			// better than the current best
			if (j == 0 || fitnesses[j] < fBest) {

				// new best fitness value
				fBest = fitnesses[j];

				// PEUDO CODE: Xbest <-- fittest individual belnoging to Pop_g
				// new best individual = jth individual
				for (int n = 0; n < problemDimension; n++)
					// The variable belonging to the jth individual in the nth dimension in the
					// population
					best[n] = population[j][n];
				FT.add(i, fBest);
			}

			i++;
		}

		// Temporary variables
		double[] x_j = new double[problemDimension];
		double[] mutant = new double[problemDimension];
		double[] x_off = new double[problemDimension];
		double f_xj = Double.NaN;
		double f_xoff = Double.NaN;
		// New population (Pop_g_+_1)
		double[][] newPop = new double[populationSize][problemDimension];

		// iterate
		while (i < maxEvaluations) {

			// PSEUDO CODE: for each x_j belonging to Pop_g
			for (int j = 0; j < populationSize && i < maxEvaluations; j++) {

				// current individual, x_j
				for (int n = 0; n < problemDimension; n++) {
					x_j[n] = population[j][n];
				}
				// and its fitness value
				f_xj = fitnesses[j];

				// mutation using DE/best/1 mutation straergy
				mutant = best1(population, best, F);

				// After mutation, crossover has to happen
				// Using Binomial Crossover
				x_off = binomialCrossover(x_j, mutant, CR);

				// Toroidal correction in case the solution lands outside the search space
				x_off = toro(x_off, bounds);

				f_xoff = problem.f(x_off);// Evaluate its fitness value
				i++;// Increase the counter after every fitness functional call

				// Replace the indiviual in the new population if the fitness value improves
				if (f_xoff < f_xj) {

					// current individual, x_j
					for (int n = 0; n < problemDimension; n++) {
						newPop[j][n] = x_off[n];
					}
					// and its fitness value
					fitnesses[j] = f_xoff;

					// Update the best solution and best fitness value if the fitness value is even
					// better than the best
					// fitness value
					if (f_xoff < fBest) {

						fBest = f_xoff;

						for (int n = 0; n < problemDimension; n++) {
							best[n] = x_off[n];
						}
						FT.add(i, fBest);// the new best fitness value is added to the fitness trend array
					}

					// Else, transfer the same individual and its fitness value to the new
					// population
				} else {
					for (int n = 0; n < problemDimension; n++) {
						newPop[j][n] = x_j[n];
					}
					fitnesses[j] = f_xj;
				}

			}
			// Make the new population equal to the current popuilation for the next
			// iteration
			population = newPop;
		}

		finalBest = best;// save the final best

		FT.add(i, fBest); // the best fitness value is added to the fitness trend array
		return FT;// return the fitness trend

	}

	// best/1 has been chosen as the mutation stratergy
	public double[] best1(double[][] population, double[] xBest, double F) {

		int problemDimension = population[0].length; // Number of dimensions of the problem
		int populationSize = population.length; // number of individuals (M) in the population

		// An array for individuals of the population
		int[] individuals = new int[populationSize];

		// transfer all individuals of the population to the array
		for (int i = 0; i < populationSize; i++)
			individuals[i] = i;

		// randomly ordered individuals of the population
		individuals = RandUtils.randomPermutation(individuals);

		// 2 randomly chosen individuals in the population,
		// with which mutaion has to happen
		int x_r1 = individuals[0];
		int x_r2 = individuals[1];

		double[] mutant = new double[problemDimension];// variable for the mutant

		// PSEUDO CODE: x_m=x_best + F(x_r1 -x_r1)
		// perform the mutation
		for (int i = 0; i < problemDimension; i++)
			mutant[i] = xBest[i] + F * (population[x_r1][i] - population[x_r2][i]);

		return mutant; // return the mutant

	}

	public static double[] best1(double[] xBest, double[] xr, double[] xs, double F) {
		int problemDimension = xr.length;
		double[] newPt = new double[problemDimension];
		for (int i = 0; i < problemDimension; i++)
			newPt[i] = xBest[i] + F * (xr[i] - xs[i]);

		return newPt;
	}

	// Binomial crossover
	public double[] binomialCrossover(double[] x_2, double[] x_1, double CR) {

		int n = x_2.length; // number of variables/dimensions in the individual
		double[] x_off = new double[n]; // offspring individual

		// PSEUDO CODE:index <-- rand int [1,n]
		// Randomly generated index
		int index = RandUtils.randomInteger(n - 1);

		for (int i = 0; i < n; i++) {
			// if CR is greater than the RandUtils.random(), uniform random number between 0
			// and 1 or its the ith index, we inherit the ith variable from the first parent
			if (RandUtils.random() < CR || i == index)
				x_off[i] = x_1[i];
			// else, we inherit the ith variable from the second parent
			else
				x_off[i] = x_2[i];
		}
		return x_off; // return the new offspring
	}
}
