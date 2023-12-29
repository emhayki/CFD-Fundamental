#include "matplotlibcpp.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

namespace plt = matplotlibcpp;

////////////////////////////////////////////////////////////
// Step 3: 1D Diffusion
////////////////////////////////////////////////////////////

int main() {
    // Simulation parameters
    const int X = 41;                    // Number of spatial points
    const int T = 20;                    // Number of time steps

    const double  v = 0.1;               // Diffusion coefficient
    const double dx = 2.0 / (X - 1);     // Spatial step size
    const double dt = 0.2 * (dx * dx) / v;    // Time step size

    // Initialize spatial grid and initial condition
    std::vector<double> x(X), u(X), un(X);

    for (int i = 0; i < X; i++) {
        x[i] = (2.0 * i) / (X - 1);
        u[i] = 1 + (x[i] >= 0.5 && x[i] <= 1);
    }

    plt::ion();
    plt::Plot plot;

    // Time-stepping loop
    for (int n = 0; n < T; n++) {
        un = u;

        for (int i = 1; i < X-1; i++) {
            u[i] = v * dt/(dx * dx) * (un[i+1] - 2*un[i] + un[i-1]) + un[i];
        }

        plot.update(x, u);
        plt::xlim(0, 2);
        plt::ylim(0.5, 2.5);
        plt::pause(0.1);
    }

    plt::show();
    return 0;
}
