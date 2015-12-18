#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
const double pi = 4*atan(1.0);

double f(double x, double t) {
    return exp(-x*x/2/t)/sqrt(2*pi*t);
}

int main() {
    ofstream file;
    char name[4][10] = {"0.01", "0.10", "1.00", "10.0"};
    double t[4] = {0.01, 0.1, 1.0, 10.0};
    for (int i = 0; i < 4; i++) {
        file.open(name[i]);
        for (int j = 0; j < 200; j++) {
            double x = -5 + j*0.005*10;
            file << x << '\t' << f(x,t[i]) << '\n';
        }
        file.close();
    }
    file.close();
}
