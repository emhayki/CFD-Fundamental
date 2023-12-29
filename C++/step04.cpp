#include "matplotlibcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace plt = matplotlibcpp;

////////////////////////////////////////////////////////////
// Step 4: 1D Burgers' Equation
////////////////////////////////////////////////////////////

double AnalyticalSolution(double t, double v, double x) {
    double phi = exp(-pow(4*t - x, 2) / (4*v*(1 + t))) 
               + exp(-pow(2*M_PI + 4*t - x, 2) / (4*v*(1 + t)));

    double dphi = (exp(-pow(2*M_PI + 4*t - x, 2) / (4*v*(1 + t))) * (4*M_PI + 8*t - 2*x)) / (4*v*(1 + t)) 
                + (exp(-pow(4*t - x, 2) / (4*v*(1 + t))) * (8*t - 2*x)) / (4*v*(1 + t));

    double u = 4 - (2 * dphi * v) / phi;

    return u;
}

int main() {
    // Simulation parameters
    const int X = 101;                         // Number of spatial points
    const int T = 100;                         // Number of time steps

    const double  v = 0.07;                    // Diffusion coefficient
    const double dx = 2. * M_PI / (X - 1);     // Spatial step size
    const double dt = dx * v;                  // Time step size

    // Initialize spatial grid and initial condition
    std::vector<double> x(X), u(X), un(X);

    for (int i = 0; i < X; i++) {
        x[i] = (2 * M_PI * i) / (X - 1);
        u[i] = AnalyticalSolution(0, v, x[i]);
    }

    plt::ion();
    plt::Plot plot;

    // Time-stepping loop
    for (int n = 0; n < T; n++) {
        un = u;

        for (int i = 1; i < X-1; i++) {
            u[i] = (v * dt / (dx * dx) * (un[i+1] - 2*un[i] + un[i-1])) - (un[i] * dt / dx * (un[i] - un[i-1])) + un[i];
        }

        // Boundary conditions
        u[0] = (v * dt / (dx * dx) * (un[1] - 2*un[0] + un[X-2])) - (un[0] * dt / dx * (un[0] - un[X-2])) + un[0];
        u[X-1] = u[0];

        plot.update(x, u);
        plt::xlim(0, 6);
        plt::ylim(0.5, 8.0);
        plt::pause(0.1);
    }

plt::show();
    return 0;
}
