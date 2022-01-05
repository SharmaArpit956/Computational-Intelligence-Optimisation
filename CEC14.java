package experiments;

import interfaces.Experiment;
import interfaces.Algorithm;
import benchmarks.CEC2014;
import algorithms.CMAES;
import algorithms.ISPO;
import algorithms.S;
import algorithms.DE;
import algorithms.jDES;

public class CEC14 extends Experiment {

	public CEC14(int probDim) throws Exception {
		// super(probDim,"cec2015allDim");
		super(probDim, 5000, "testCEC14");
		setNrRuns(30);

		Algorithm a;// ///< A generic optimiser.

		a = new jDES(); // make an instance of jDES algorithm 
		a.setParameter("p0", 0.1);// F_lowerBound
		a.setParameter("p1", 1.0);// F_upperBound
		a.setParameter("p2", 0.1);// tau1
		a.setParameter("p3", 0.1);// tau2
		a.setParameter("p4", 150.0);// maximum local iterations
		a.setParameter("p5", 0.15);// alpha
		add(a);// add it to the list

		a = new DE(); // make an instance of jDES algorithm  
		a.setParameter("p0", 0.5);// F
		a.setParameter("p1", 0.8);// CR
		add(a);// add it to the list

		a = new S(); // make an instance of S algorithm defined in S.java file
		add(a);// add it to the list

		for (int i = 1; i <= 30; i++)
			add(new CEC2014(probDim, i));

	}
}
