// Fokker-Planck equation approach to VMC

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "rng.h"

using namespace std;

int N;                          // number of walkers
double *x;                      // positions of walkers
double alpha;                   // variational parameter
double tStep;                   // time step
int seed = -987654321;          // for random ran2 and gasdev

double eSum;                    // accumulator to find energy
double eSqdSum;                 // accumulator to find fluctuations in E

void initialize() {

    x = new double [N];
    for (int i = 0; i < N; i++)
        x[i] = qadran() - 0.5;
    tStep = 0.1;
}

void zeroAccumulators() {
    eSum = eSqdSum = 0;
}

double eLocal(double x) {

    // compute the local energy
    return alpha + x * x * (0.5 - 2 * alpha * alpha);
}

int nAccept;                   // accumulator for number of accepted steps

void MetropolisStep(int n) {

    // make a trial move
    double x = ::x[n];              // :: chooses the global x
    double Fx = - 4 * alpha * x;
    double y = x + gasdev(seed) * sqrt(tStep) + Fx * tStep / 2;
    // compute ratio for Metropolis test
    double rhoRatio = exp( - 2 * alpha * (y * y - x * x));
    double oldExp = y - x - Fx * tStep / 2;
    double Fy = - 4 * alpha * y;
    double newExp = x - y - Fy * tStep / 2;
    double GRatio = exp( -(newExp * newExp - oldExp * oldExp) / (2 * tStep));
    double w = rhoRatio * GRatio;

    // Metropolis test
    if (w > ran2(seed)) {
        ::x[n] = x = y;
        ++nAccept;
    }

    // accumulate energy and wave function
    double e = eLocal(x);
    eSum += e;
    eSqdSum += e * e;
}

void oneMonteCarloStep() {

    // perform N Metropolis steps
    for (int n = 0; n < N; n++) {
        MetropolisStep(n);
    }
}

int main() {

    cout << " Fokker-Planck approach to VMC: Harmonic Oscillator\n"
         << " --------------------------------------------------\n"
         << " Enter number of walkers: ";
    cin >> N;
    cout << " Enter variational parameter alpha: ";
    cin >> alpha;
    cout << " Enter number of Monte Carlo steps: ";
    int MCSteps;
    cin >> MCSteps;

    initialize();

    // perform 20% of MCSteps as thermalization steps
    // and adjust time step size so acceptance ratio ~90%
    int thermSteps = int(0.2 * MCSteps);
    int adjustInterval = int(0.1 * thermSteps) + 1;
    nAccept = 0;
    cout << " Performing " << thermSteps << " thermalization steps ..." 
         << flush;
    for (int i = 0; i < thermSteps; i++) {
        oneMonteCarloStep();
        if ((i+1) % adjustInterval == 0) {
            tStep *= nAccept / (0.9 * N * adjustInterval);
            nAccept = 0;
        }
    }
    cout << "\n Adjusted time step size = " << tStep << endl;

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
}
