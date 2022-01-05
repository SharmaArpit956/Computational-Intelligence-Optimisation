package algorithms;

 import static utils.algorithms.Misc.generateRandomSolution;
import interfaces.Algorithm;
import interfaces.Problem;
import utils.RunAndStore.FTrend;
import static utils.algorithms.Misc.*;
import utils.random.RandUtils;

 public class jDES extends Algorithm {

	@Override
	public FTrend execute(Problem problem, int maxEvaluations) throws Exception {

		// Scaling Factor for the difference vector for the mutatioon
		double F_lowerBound = this.getParameter("p0");// 0.1
		double F_upperBound = this.getParameter("p1");// 1
		double tau1 = this.getParameter("p2"); // 0.1
		double tau2 = this.getParameter("p3"); // 0.1
		int maximumLocalIterations = getParameter("p4").intValue(); // maximum local iterations

		// alpha for calculating the exploratory radius of the S algorithm, local
		// searcher
		double alpha = getParameter("p5").doubleValue();

		FTrend FT = new FTrend();
		int problemDimension = problem.getDimension(); // Number of dimensions of the problem

		int populationSize = problemDimension * 10; // variable 'M' in the Pseudo code

		double[][] bounds = problem.getBounds(); // An array containing boundaries for each dimension

		// An 2d Array containing population in every dimension of the problem
		double[][] population = new double[populationSize][problemDimension];

	 
		double[] F = new double[populationSize];
		double[] CR = new double[populationSize];
	 
		for (int j = 0; j < populationSize; j++) {
			F[j] = F_lowerBound + RandUtils.random() * (F_upperBound - F_lowerBound);
			CR[j] = RandUtils.random();
		}

		// Initialize the varible for the fitness values of each invidual of the
		// population of the corresponding
		// dimension
		double[] fitnesses = new double[populationSize];

		// initialize the fittest individual
		double[] best = new double[problemDimension];

		// initialize the best fitness value
		double fBest = Double.NaN;

		int i = 0;

		// evaluate initial population
		for (int j = 0; j < populationSize; j++) {

			// PSEUDO CODE: Pop_g <-- randomly sample M n-dimentional individuals within D
			// Generate random individuals for the entire population within the search
			// space
			double[] tmp = generateRandomSolution(bounds, problemDimension);
			for (int n = 0; n < problemDimension; n++)
				population[j][n] = tmp[n];

			// The fitness value of each invidual of the population
			fitnesses[j] = problem.f(population[j]);

			// if its the first individual or the fitness value of the jth individual is
			// better than the current best
			if (j == 0 || fitnesses[j] < fBest) {

				// new best fitness value
				fBest = fitnesses[j];

				// PEUDO CODE: Xbest <-- fittest individual belnoging to Pop_g
				// new best individual
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

				// TODO
				// update F
				if (RandUtils.random() < tau1)
					F[j] = F_lowerBound + F_upperBound * RandUtils.random();

				// mutation using DES/rand/2 mutation straergy
				mutant = rand2(population, F[j]);

				// TODO
				// update CR
				if (RandUtils.random() < tau2)
					CR[j] = RandUtils.random();

				// After mutation, crossover has to happen
				// Using Binomial Crossover
				x_off = binomialCrossover(x_j, mutant, CR[j]);

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

						// Perform the local search in the vicinity using S as a local searcher
						// as there are high chances of the the best solution to lie near the current
						// best
						/////////////////////////////////////////////////////////
						//////////////////////////// S STARTS HERE /////////////

						int i_local = 0; // The counter for keeeping track of the local computational budget consumed
						double[] x_short = new double[problemDimension];// solution for the S algorithm

						// the best solution till now is the first solution as we donot have any other
						// solution yet
						for (int n = 0; n < problemDimension; n++) {
							best[n] = x_off[n];
							x_short[n] = x_off[n];
						}

						// Similarly, the best fitness value till now is the first finess value as we
						// donot any other finess value yet
						fBest = f_xoff;

						// temp variables
						double[] radius_local = new double[problemDimension]; // exploratory radius

						// makes exploratory radius equal to α(Upper limit -
						// lower limit) for every component of the solution
						for (int p_local = 0; p_local < problemDimension; p_local++) {
							radius_local[p_local] = alpha * (bounds[p_local][1] - bounds[p_local][0]);
						}

						// Flag for determining if the fitness ever improved or not. As no fitness value
						// improvement has been made till now, it is initialised as true
						boolean neverImproved_local = true;

						// S Algorithm here
						while (i_local < maximumLocalIterations && i < maxEvaluations) { // keep the loop running till
																							// the budget expires

							neverImproved_local = true; // re-initialize the flag to true for the next iteration

							// for every axis of the variable and till the budget expires
							for (int j_local = 0; j_local < problemDimension && i < maxEvaluations; j_local++) {

								// perturb the jth component (starting from 0) towards negative jth axis
								// treating other components as constants (like partial derivative)
								x_short[j_local] = best[j_local] - radius_local[j_local];

								// After every perturbation of the solution,a correction method is needed as as
								// there are chances of the component of having landed outside the search space
								x_short = toro(x_short, bounds);

								// fitness value calculaion acoording to the problem being handled
								f_xoff = problem.f(x_short);

								i++; // Increase the counter after every fitness functional call

								// Fitness Improved (Found a lower point on the curve for the jth dimension)
								if (f_xoff <= fBest) {

									// Very important. this means that improvement in fitness value happened.
									// Therefore, we cannot say that improvement never happened
									neverImproved_local = false;

									best[j_local] = x_short[j_local]; // As the improvenment happened, the best
																		// solution is updated
									fBest = f_xoff; // Similarly, the best fitness value is also updated

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
									x_short[j_local] = best[j_local] + (radius_local[j_local] / 2);

									// After every perturbation of the solution, correction method is needed as as
									// there are chances of x landing outside the search space in both cases
									x_short = toro(x_short, bounds);

									f_xoff = problem.f(x_short); // fitness functional call
									i++; // Increase the counter after every fitness functional call

									// Fitness Improved (Found a lower point on the curve for the jth dimension)
									if (f_xoff <= fBest) {

										// Very important. this means that improvement in fitness value happened.
										// Therefore, we cannot say improvement never happened
										neverImproved_local = false;

										best[j_local] = x_short[j_local]; // As the improvenment happened, the
																			// best solution is updated
										fBest = f_xoff; // Similarly, the best fitness value is also updated

										// the new best fitness value is added to the fitness trend array
										FT.add(i, fBest);

									} else { // if the fitness value doesnt improve
										x_short[j_local] = best[j_local]; // store the best solution as the
																			// current solution

										// store the best fitness functions as the current fitness function
										f_xoff = fBest;
									}
								}
							}

							// If the fitness function never improved throughout the for loop
							if (neverImproved_local == true) {

								for (int d_local = 0; d_local < problemDimension; d_local++) {

									// half the radius to do the perturbations as it means there
									// are more chances of the minima lying near the current point rather
									// than far from it(which we have checked and not found at this point)
									radius_local[d_local] = radius_local[d_local] / 2;
								}
							}
							i_local++; // increment local iterations counter after every cycle
						}

						// make x_off equal to x_short after the local search
						for (int n = 0; n < problemDimension; n++) {
							x_off[n] = x_short[n];
						}
						/////////////////////////////////////////////////////////
						//////////////////////////// S ENDS HERE/////////////////

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

	// rand2 has been chosen as the mutation stratergy
	public double[] rand2(double[][] population, double F) {

		int problemDimension = population[0].length; // Number of dimensions of the problem
		int populationSize = population.length; // number of individuals (M) in the population

		// An array for individuals of the population
		int[] individuals = new int[populationSize];

		// transfer all individuals of the population to the array
		for (int i = 0; i < populationSize; i++)
			individuals[i] = i;

		// randomly ordered individuals of the population
		individuals = RandUtils.randomPermutation(individuals);

		// five randomly chosen individuals in the population,
		// with which mutaion has to happen
		int x_r1 = individuals[0];
		int x_r2 = individuals[1];
		int x_r3 = individuals[2];
		int x_r4 = individuals[3];
		int x_r5 = individuals[4];

		double[] mutant = new double[problemDimension];// variable for the mutant

		// PSEUDO CODE: x_m=x_r1 + F(x_r2 -x_r3) + F(x_r4 -x_r5)
		// perform the mutation
		for (int i = 0; i < problemDimension; i++)
			mutant[i] = population[x_r1][i] + F * (population[x_r2][i] - population[x_r3][i])
					+ F * (population[x_r4][i] - population[x_r5][i]);

		return mutant; // return the mutant
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
