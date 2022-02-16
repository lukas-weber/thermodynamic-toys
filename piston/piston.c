#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

double frand() {
	return random()/(double)RAND_MAX;
}

// taken from https://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/
double randn() {
	double U1, U2, W, mult;
	static double X1, X2;
	static int call = 0;

	if (call == 1) {
		call = !call;
		return X2;
	}

	do {
		U1 = -1 + frand()*2;
		U2 = -1 + frand()*2;
		W = pow(U1, 2) + pow(U2, 2);
	} while (W >= 1 || W == 0);

	mult = sqrt((-2*log(W)) / W);
	X1 = U1*mult;
	X2 = U2*mult;

	call = !call;

	return X1;
}


enum {
	EDGE_NONE, // should not happen
	EDGE_LEFT,
	EDGE_RIGHT
};

void check(int Npart, double *xs, double X) {
	double xmin = INFINITY;
	double xmax = -INFINITY;
	for(int i = 0; i < Npart; i++) {
		if(xs[i] > xmax)
			xmax = xs[i];
		if(xs[i] < xmin)
			xmin = xs[i];
	}

	printf("0 < %g < %g < %g\n", xmin, xmax, X);
	assert(0 < xmin);
	assert(xmax <= X);
}

void simulate_piston(double *outPositions, int timeSteps, double deltat, int Npart, double mpart, double T, double X0, double gravityg) {
	double L0 = Npart*T/gravityg;

	// particle coordinates
	double *xs = malloc(sizeof(double)*Npart);
	double *vs = malloc(sizeof(double)*Npart);

	// piston coordinates
	double X = L0*X0;
	double V = 0;

	double time = 0;
	int step = 1;
	assert(timeSteps > 0);
	outPositions[0] = X/L0;

	for(int i = 0; i < Npart; i++) {
		xs[i] = L0*frand();
		vs[i] = sqrt(T/mpart)*randn(); 
	}

	while(step < timeSteps) {
		double minTime = INFINITY;
		int collidingParticle = -1;
		int edge = EDGE_NONE;
		// find next collision
		for(int i = 0; i < Npart; i++) {
			double timeToLeft = (0-xs[i])/(vs[i]);
			if(timeToLeft > 0 && timeToLeft < minTime) {
				minTime = timeToLeft;
				edge = EDGE_LEFT;
				collidingParticle = i;
			}

			// solve quadratic equation
			// x + v*t = X + V*t - g*t^2/2
			double p = (V-vs[i])/gravityg;
			double q = p*p + 2*(X-xs[i])/gravityg;

			// Solution exists
			if(q > 0) { 
				q = sqrt(q);
				double timeToRight = p - q;
				if(timeToRight < 0)
					timeToRight += 2*q;

				if(timeToRight > 0 && timeToRight < minTime) {
					minTime = timeToRight;
					edge = EDGE_RIGHT;
					collidingParticle = i;
				}
			}
		}

		while(step < timeSteps && (step + 1)*deltat < time + minTime) {
			double timeToStep = step*deltat-time;
			outPositions[step] = (X + V * timeToStep - gravityg*timeToStep*timeToStep/2)/L0;
			step++;
		}

		// move to collision point
		for(int i = 0; i < Npart; i++) {
			xs[i] += minTime*vs[i];
		}
		
		X += V*minTime - gravityg*minTime*minTime/2;
		time += minTime;
		V -= gravityg*minTime;

		assert(edge != EDGE_NONE);

		//printf("colliding: %d, edge: %d,  minTime: %g\n", collidingParticle, edge, minTime);
		//printf("(%g %g) (%g %g)\n", xs[collidingParticle], X, vs[collidingParticle], V);

		if(edge == EDGE_LEFT) {
			vs[collidingParticle] *= -1;
		} else {
			// elastic collision formula
			double vnew = (vs[collidingParticle]*(mpart-1)+2*V)/(mpart+1);
			double Vnew = (V*(1-mpart)+2*mpart*vs[collidingParticle])/(mpart+1);
			V = Vnew;
			vs[collidingParticle] = vnew;
		}


		// regularization to compensate for rounding errors
		// 
		// propagate a bit further so we don’t repeat the same collision
		double epsilon = 1e-8;
		double xmax = 0;
		for(int i = 0; i < Npart; i++) {
			xs[i] += epsilon*vs[i];
			if(xs[i] > xmax)
				xmax = xs[i];
		}
		X += V*epsilon-gravityg*epsilon*epsilon/2;
		V -= gravityg*epsilon;
		// in the end everything has to be inside of the piston.
		if(X < xmax)
			X = xmax+epsilon;

		//printf("-> (%g %g) (%g %g)\n", xs[collidingParticle], X, vs[collidingParticle], V);
		//check(Npart, xs, X);
	}
}
			

