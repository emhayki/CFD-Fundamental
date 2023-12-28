#include "matplotlibcpp.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

namespace plt = matplotlibcpp;

////////////////////////////////////////////////////////////
// Step 2: 1D Nonlinear Convection
////////////////////////////////////////////////////////////

int main() {
    // Simulation parameters
    const int X = 41;                    // Number of spatial points
    const int T = 30;                    // Number of time steps

    const double dx = 2.0 / (X - 1);     // Spatial step size
    const double dt = 0.2 * dx;          // Time step size

    // Initialize spatial grid and initial condition
    std::vector<double> x(X), u(X), un(X);

    for (int i = 0; i < X; i++) {
        x[i] = (5.0 * i) / (X - 1);
        u[i] = (x[i] >= 0.5 && x[i] <= 1) ? 2 : 1;
    }

    plt::ion();
    plt::Plot plot;

    // Time-stepping loop
    for (int n = 0; n < T; n++) {
        un = u;

        for (int i = 1; i < X; i++) {
            u[i] = un[i] - un[i] * (un[i] - un[i-1]) * dt / dx;
        }

        plot.update(x, u);
        plt::xlim(0, 2);
        plt::ylim(0.5, 2.5);
        plt::pause(0.1);
    }

    plt::show();
    return 0;
}
