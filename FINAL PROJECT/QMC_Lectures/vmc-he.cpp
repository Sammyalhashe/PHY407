// Variational Monte Carlo for the Helium Atom

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "rng.h"

using namespace std;

const int NDIM = 3;        // dimensionality of space
const int NELE = 2;        // number of electrons
int N;                     // number of walkers
double (*r)[NELE][NDIM];   // walker coordinates in 6-D configuration space

double alpha;              // Pade-Jastrow variational parmeter
double delta;              // trial step size

void initialize() {
    r = new double [N][NELE][NDIM]; 
    for (int n = 0; n < N; n++)
    for (int e = 0; e < NELE; e++)
    for (int d = 0; d < NDIM; d++)
        r[n][e][d] = qadran() - 0.5;
    delta = 1;
}

double eSum;
double eSqdSum;

void zeroAccumulators() {
    eSum = eSqdSum = 0;
}

double Psi(double *rElectron1, double *rElectron2) {

    // value of trial wave function for walker n
    double r1 = 0, r2 = 0, r12 = 0;
    for (int d = 0; d < 3; d++) {
        r1 += rElectron1[d] * rElectron1[d];
        r2 += rElectron2[d] * rElectron2[d];
        r12 += (rElectron1[d] - rElectron2[d]) 
               * (rElectron1[d] - rElectron2[d]);
    }
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    r12 = sqrt(r12);
    double Psi = - 2*r1 - 2*r2 + r12 / (2 * (1 + alpha*r12));
    return exp(Psi);

}

double eLocal(double *rElectron1, double *rElectron2) {

    // value of trial wave function for walker n
    double r1 = 0, r2 = 0, r12 = 0;
    for (int d = 0; d < 3; d++) {
        r1 += rElectron1[d] * rElectron1[d];
        r2 += rElectron2[d] * rElectron2[d];
        r12 += (rElectron1[d] - rElectron2[d]) *
               (rElectron1[d] - rElectron2[d]);
    }
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    r12 = sqrt(r12);
    double dotProd = 0;
    for (int d = 0; d < 3; d++) {
        dotProd += (rElectron1[d] - rElectron2[d]) / r12 *
                   (rElectron1[d] / r1 - rElectron2[d] / r2);
    }
    double denom = 1 / (1 + alpha * r12);
    double denom2 = denom * denom;
    double denom3 = denom2 * denom;
    double denom4 = denom2 * denom2;
    double e = - 4 + alpha * (denom + denom2 + denom3) 
               - denom4 / 4 + dotProd * denom2;
    return e;
}

int nAccept;

void MetropolisStep(int walker) {

    // make a trial move of each electron
    double rElectron1[3], rElectron2[3], rTrial1[3], rTrial2[3];
    for (int d = 0; d < 3; d++) {
        rElectron1[d] = r[walker][0][d];
        rTrial1[d] = rElectron1[d] + delta * (2 * qadran() - 1);
        rElectron2[d] = r[walker][1][d];
        rTrial2[d] = rElectron2[d] + delta * (2 * qadran() - 1);
    }

    // Metropolis test
    double w = Psi(rTrial1, rTrial2) / Psi(rElectron1, rElectron2);
    if (qadran() < w * w) {
        for (int d = 0; d < 3; d++) {
            r[walker][0][d] = rElectron1[d] = rTrial1[d];
            r[walker][1][d] = rElectron2[d] = rTrial2[d];
        }
        ++nAccept;
    }

    // accumulate local energy
    double e = eLocal(rElectron1, rElectron2);
    eSum += e;
    eSqdSum += e * e;
}

void oneMonteCarloStep() {

    // do Metropolis step for each walker
    for (int n = 0; n < N; n++)
        MetropolisStep(n);
}

int main() {

    cout << " Variational Monte Carlo for Helium Atom\n"
         << " ---------------------------------------\n";
    cout << " Enter number of walkers: ";
    cin >> N;
    cout << " Enter parameter Pade-Jastrow parameter alpha: ";
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
    cout << "\n Adjusted step size delta = " << delta << endl;

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
