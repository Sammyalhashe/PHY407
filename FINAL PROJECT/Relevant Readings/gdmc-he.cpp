// Guide Function Diffusion Monte Carlo for Helium Atom

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "rng.h"

using namespace std;

int seed = -987654321;            // seed for ran2 and gasdev

int N;                            // current number of walkers
int N_T;                          // desired target number of walkers
const int DIM = 6;                // dimension of R = (r1, r2)
double **R;                       // walker positions in 6-D space
bool *alive;                      // is this walker alive?

void ensureCapacity(int index) {

    static int maxN = 0;          // remember size of the arrays

    if (index < maxN)
        return;                   // additional storage not needed

    int oldN = maxN;              // remember old capacity to copy values
    if (maxN > 0)
        maxN *= 2;                // double capacity
    else
        maxN = 1;
    if (index > maxN - 1)         // if this is not enough
        maxN = index + 1;         // increase to make it enough

    // allocate new storage
    double **newR = new double* [maxN];
    bool *newAlive = new bool [maxN];
    for (int n = 0; n < maxN; n++) {
        newR[n] = new double [DIM];
        if (n < oldN) {
            for (int d = 0; d < DIM; d++)
                newR[n][d] = R[n][d];
            newAlive[n] = alive[n];
            delete [] R[n];       // release old memory
        }
    }

    // delete old storage and point to new
    delete [] R;
    R = newR;
    delete [] alive;
    alive = newAlive;
}

double ESum, ESqdSum;             // accumulators for observables

void zeroAccumulators() {
    ESum = ESqdSum = 0;
}

double dt;                        // time step Delta_t set by user
double E_T;                       // target energy

void initialize() {

    // create target number of walkers
    N = N_T;
    for (int n = 0; n < N; n++) {
        ensureCapacity(n);
        for (int d = 0; d < DIM; d++)
            R[n][d] = ran2(seed) - 0.5;
        alive[n] = true;
    }

    // set target energy close to VMC result
    E_T = -2.85;
}

void findSeparations(double *R, double& r1, double& r2, double&r12) {

    // find electron-nucleus and electron-electron separations
    r1 = r2 = r12 = 0;
    for (int e1 = 0; e1 < 3; e1++) {
        int e2 = e1 + 3;           // second electron indices
        r1 += R[e1] * R[e1];
        r2 += R[e2] * R[e2];
        r12 += (R[e1] - R[e2]) * (R[e1] - R[e2]);
    }
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    r12 = sqrt(r12);
}

double alpha = 0.15;              // Pade-Jastrow wave function parameter

double Psi_T(double *R) {

    // value of guide function
    double r1, r2, r12;
    findSeparations(R, r1, r2, r12);
    double Psi_T = - 2 * r1 - 2 * r2 + r12 / (2 * (1 + alpha * r12));
    return exp(Psi_T);
}

double E_L(double *R) {

    // value of local energy for guide function
    double r1, r2, r12;
    findSeparations(R, r1, r2, r12);
    double dotProd = 0;
    for (int e1 = 0; e1 < 3; e1++) {
        int e2 = e1 + 3;          // second electron indices
        dotProd += (R[e1] - R[e2]) / r12 * (R[e1] / r1 - R[e2] / r2);
    }
    double denom = 1 / (1 + alpha * r12);
    double denom2 = denom * denom;
    double denom3 = denom2 * denom;
    double denom4 = denom2 * denom2;
    double E_L = - 4 + alpha * (denom + denom2 + denom3) 
               - denom4 / 4 + dotProd * denom2;
    return E_L;
}

void findForce(double *R, double *F) {

    // find Fokker-Planck forces
    double r1, r2, r12;
    findSeparations(R, r1, r2, r12);
    for (int d = 0; d < DIM; d++)
        F[d] = 0;
    double denom2 = 1 / (1 + alpha * r12);
    denom2 *= denom2;
    for (int e1 = 0; e1 < 3; e1++) {
        int e2 = e1 + 3;
        F[e1] += - 4 * R[e1] / r1 + denom2 * (R[e1] - R[e2]) / r12;
        F[e2] += - 4 * R[e2] / r2 + denom2 * (R[e2] - R[e1]) / r12;
    }
}

double G(double *RPrime, double *R) {

    // value of Fokker-Planck Green's function exponential
    double F[DIM];
    findForce(R, F);
    double G = 0;
    for (int d = 0; d < DIM; d++) {
        double dR = RPrime[d] - R[d] - F[d] * dt / 2;
        G += dR * dR;
    }
    return exp(- G / (2 * dt));
}

int nTrials;                    // number of Metropolis tests
int nAccept;                    // number of acceptances

void oneMonteCarloStep(int n) {

    // define position variables for this walker
    double R[DIM];
    for (int d = 0; d < DIM; d++) 
        R[d] = ::R[n][d];
    // trial shift walker to new position
    double F[DIM], RPrime[DIM];
    findForce(R, F);
    for (int d = 0; d < DIM; d++) 
        RPrime[d] = R[d] + F[d] * dt / 2 + gasdev(seed) * sqrt(dt);
    // Metropolis acceptance test
    double w = Psi_T(RPrime) / Psi_T(R);
    w *= w * G(R, RPrime) / G(RPrime, R);
    ++nTrials;
    if (w > ran2(seed)) {
        for (int d = 0; d < DIM; d++) 
            ::R[n][d] = RPrime[d];
        ++nAccept;
    } else {
        return;                 // don't do branching step below
    }
    // branching step
    double q = exp(- dt * (E_L(RPrime) - E_T));  

    int survivors = int(q);
    if (q - survivors > ran2(seed))
        ++survivors;

    // append survivors-1 copies of the walker to the end of the array
    for (int i = 0; i < survivors - 1; i++) {
        ensureCapacity(N);
        for (int d = 0; d < DIM; d++)
            ::R[N][d] = RPrime[d];
        alive[N] = true;
        ++N;
    }

    // if survivors is zero, then kill the walker
    if (survivors == 0)
        alive[n] = false;
}

void oneTimeStep() {

    // DMC step for each walker
    int N_0 = N;
    for (int n = 0; n < N_0; n++)
        oneMonteCarloStep(n);

    // adjust E_T
    E_T += log(N_T / double(N)) / 10;

    // remove all dead walkers from the arrays
    int newN = 0;
    for (int n = 0; n < N; n++)
    if (alive[n]) {
        if (n != newN) {
            for (int d = 0; d < DIM; d++)
                R[newN][d] = R[n][d];
            alive[newN] = true;
        }
        ++newN;
    }
    N = newN;

    // measure energy, wave function
    ESum += E_T;
    ESqdSum += E_T * E_T;
}

int main() {

    cout << " Guide Function Diffusion Monte Carlo for Helium Atom\n"
         << " ----------------------------------------------------\n";
    cout << " Enter desired target number of walkers: ";
    cin >> N_T;
    cout << " Enter time step dt: ";
    cin >> dt;
    cout << " Enter total number of time steps: ";
    int timeSteps;
    cin >> timeSteps;

    initialize();

    // do 20% of timeSteps as thermalization steps
    int thermSteps = int(0.2 * timeSteps);
    int adjustInterval = int(0.1 * thermSteps) + 1;
    nTrials = nAccept = 0;
    cout << " Performing " << thermSteps << " thermalization steps ..." 
         << flush;
    for (int i = 0; i < thermSteps; i++) {
        oneTimeStep();
        if ((i+1) % adjustInterval == 0) {
            dt *= nAccept / (0.9 * nTrials);
            nTrials = nAccept = 0;
        }
    }
    cout << "\n Adjusted time step size = " << dt << endl;

    // production steps
    zeroAccumulators();
    for (int i = 0; i < timeSteps; i++) {
        oneTimeStep();
    }

    // compute averages
    double EAve = ESum / timeSteps;
    double EVar = ESqdSum / timeSteps - EAve * EAve;
    cout << " <E> = " << EAve << " +/- " << sqrt(EVar / timeSteps) << endl;
    cout << " <E^2> - <E>^2 = " << EVar << endl;
}
