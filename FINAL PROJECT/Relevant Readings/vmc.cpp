// Variational Monte Carlo for the harmonic oscillator

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "rng.h"
using namespace std;

int seed = -123456789;         // for random number generator

int N;                         // number of walkers
double *x;                     // walker positions
double delta;                  // step size

double eSum;                   // accumulator to find energy
double eSqdSum;                // accumulator to find fluctuations in E
double xMin = -10;             // minimum x for histogramming psi^2(x)
double xMax = +10;             // maximum x for histogramming psi^2(x)
double dx = 0.1;               // psi^2(x) histogram step size
double *psiSqd;                // psi^2(x) histogram
int nPsiSqd;                   // size of array

void zeroAccumulators() {
    eSum = eSqdSum = 0;
    for (int i = 0; i < nPsiSqd; i++)
        psiSqd[i] = 0;
}

void initialize() {

    x = new double [N];
    for (int i = 0; i < N; i++)
        x[i] = ran2(seed) - 0.5;
    delta = 1;

    nPsiSqd = int((xMax - xMin) / dx);
    psiSqd = new double [nPsiSqd];

    zeroAccumulators();
}

double alpha;                  // trial function is exp(-alpha*x^2)

double p(double xTrial, double x) {

    // compute the ratio of rho(xTrial) / rho(x)
    return exp(- 2 * alpha * (xTrial*xTrial - x*x));
}

double eLocal(double x) {

    // compute the local energy
    return alpha + x * x * (0.5 - 2 * alpha * alpha);
}

int nAccept;                   // accumulator for number of accepted steps

void MetropolisStep() {

    // chose a walker at random
    int n = int(ran2(seed) * N);

    // make a trial move
    double xTrial = x[n] + delta * gasdev(seed);

    // Metropolis test
    if (p(xTrial, x[n]) > ran2(seed)) {
        x[n] = xTrial;
        ++nAccept;
    }

    // accumulate energy and wave function
    double e = eLocal(x[n]);
    eSum += e;
    eSqdSum += e * e;
    int i = int((x[n] - xMin) / dx);
    if (i >= 0 && i < nPsiSqd)
        psiSqd[i] += 1;
}

void oneMonteCarloStep() {

    // perform N Metropolis steps
    for (int i = 0; i < N; i++) {
        MetropolisStep();
    }
}

int main() {

    cout << " Variational Monte Carlo for Harmonic Oscillator\n"
         << " -----------------------------------------------\n";
    cout << " Enter number of walkers: ";
    cin >> N;
    cout << " Enter parameter alpha: ";
    cin >> alpha;
    cout << " Enter number of Monte Carlo steps: ";
    int MCSteps;
    cin >> MCSteps;

    initialize();

    // perform 20% of MCSteps as thermalization steps
    // and adjust step size so acceptance ratio ~50%
    int thermSteps = int(0.2 * MCSteps);
    int adjustInterval = int(0.1 * thermSteps) + 1;
    nAccept = 0;
    cout << " Performing " << thermSteps << " thermalization steps ..." 
         << flush;
    for (int i = 0; i < thermSteps; i++) {
        oneMonteCarloStep();
        if ((i+1) % adjustInterval == 0) {
            delta *= nAccept / (0.5 * N * adjustInterval);
            nAccept = 0;
        }
    }
    cout << "\n Adjusted Gaussian step size = " << delta << endl;

    // production steps
    zeroAccumulators();
    nAccept = 0;
    cout << " Performing " << MCSteps << " production steps ..." << flush;
    for (int i = 0; i < MCSteps; i++)
        oneMonteCarloStep();

    // compute and print energy
    double eAve = eSum / double(N) / MCSteps;
    double eVar = eSqdSum / double(N) / MCSteps - eAve * eAve;
    double error = sqrt(eVar) / sqrt(double(N) * MCSteps);
    cout << "\n <Energy> = " << eAve << " +/- " << error
         << "\n Variance = " << eVar << endl;

    // write wave function squared in file
    ofstream file("psiSqd.data");
    double psiNorm = 0;
    for (int i = 0; i < nPsiSqd; i++)
        psiNorm += psiSqd[i] * dx;
    for (int i = 0; i < nPsiSqd; i++) {
        double x = xMin + i * dx;
        file << x << '\t' << psiSqd[i] / psiNorm << '\n';
    }
    file.close();
    cout << " Probability density written in file psiSqd.data" << endl;
}
