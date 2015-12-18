// Variational Monte Carlo for the Helium Atom with Z and beta

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "rng.h"

using namespace std;

const int NDIM = 3;        // dimensionality of space
const int NELE = 2;        // number of electrons
int N;                     // number of walkers
double (*r)[NELE][NDIM];   // walker coordinates in 6-D configuration space

double alpha;              // Pade-Jastrow variational parmeter
double delta;              // trial step size
double Z = 2;              // effective nuclear charge parameter
double beta = 0.5;         // effective electron-electron coupling parameter

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
    double Psi = - Z*r1 - Z*r2 + beta * r12 / (1 + alpha*r12);
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
    double e = - Z * Z + (Z - 2) *(1 / r1 + 1 / r2)
               + 1 / r12 * (1 - 2 * beta * denom2) + 2 * alpha * beta * denom3
               - beta * beta * denom4 + Z * beta * dotProd * denom2;
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

double eAve;             // average energy
double eVar;             // variance in the energy
int MCSteps = 10000;     // number of Monte Carlo steps per walker

void runMonteCarlo() {

    // perform 20% of MCSteps as thermalization steps
    // and adjust step size so acceptance ratio ~50%
    int thermSteps = int(0.2 * MCSteps);
    int adjustInterval = int(0.1 * thermSteps) + 1;
    nAccept = 0;
    for (int i = 0; i < thermSteps; i++) {
        oneMonteCarloStep();
        if ((i+1) % adjustInterval == 0) {
            delta *= nAccept / (0.5 * N * adjustInterval);
            nAccept = 0;
        }
    }

    // production steps
    zeroAccumulators();
    nAccept = 0;
    for (int i = 0; i < MCSteps; i++)
        oneMonteCarloStep();
    eAve = eSum / double(N) / MCSteps;
    eVar = eSqdSum / double(N) / MCSteps - eAve * eAve;
}

int main() {

    cout << " Variational Monte Carlo for Helium Atom\n"
         << " ---------------------------------------\n";
    N = 300;
    cout << " Number of walkers = " << N << endl;
    cout << " Number of MCSteps = " << MCSteps << endl;
    initialize();

    // Vary Z holding beta and alpha fixed
    ofstream file("Z.data");
    beta = 0.5;
    alpha = 0.1;
    Z = 0.5;
    cout << " Varying Z with beta = 0.5 and alpha = 0.1" << endl;
    while (Z < 3.6) {
        runMonteCarlo();
        file << Z << '\t' << eAve << '\t' << eVar << '\n';
        cout << " Z = " << Z << "\t<E> = " << eAve 
             << "\tVariance = " << eVar << endl;
        Z += 0.25;
    }
    file.close();

    // Vary beta holding Z and alpha fixed
    file.open("beta.data");
    Z = 2;
    beta = 0;
    cout << " Varying beta with Z = 2 and alpha = 0.1" << endl;
    while (beta < 1.1) {
        runMonteCarlo();
        file << beta << '\t' << eAve << '\t' << eVar << '\n';
        cout << " beta = " << beta << "\t<E> = " << eAve 
             << "\tVariance = " << eVar << endl;
        beta += 0.1;
    }
    file.close();

    // Vary alpha holding Z and beta fixed
    file.open("alpha.data");
    beta = 0.5;
    Z = 2;
    alpha = 0;
    cout << " Varying alpha with Z = 2 and beta = 0.5" << endl;
    while (alpha < 0.6) {
        runMonteCarlo();
        file << alpha << '\t' << eAve << '\t' << eVar << '\n';
        cout << " alpha = " << alpha << "\t<E> = " << eAve 
             << "\tVariance = " << eVar << endl;
        alpha += 0.05;
    }
    file.close();
}
