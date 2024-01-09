#include "matplotlibcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace plt = matplotlibcpp;

////////////////////////////////////////////////////////////
// Step 10: 2D Poisson Equation
////////////////////////////////////////////////////////////

int main() {
    // Define simulation parameters
    const int X = 50;
    const int Y = 50;
    const int T = 100;

    const double xmin = 0, xmax = 2;
    const double ymin = 0, ymax = 1;

    const double dx = (xmax - xmin) / (X - 1);
    const double dy = (ymax - ymin) / (Y - 1);

    // Create spatial grids
    std::vector<double> x(X), y(Y);

    for (int i = 0; i < X; i++)
        x[i] = xmin + i * dx;

    for (int i = 0; i < Y; i++)
        y[i] = ymin + i * dy;

    std::vector<std::vector<double>> nX(X, std::vector<double>(Y)), nY(X, std::vector<double>(Y));
    for (int i = 0; i < X; ++i) {
        for (int j = 0; j < Y; ++j) {
            nX[i][j] = x[i]; 
            nY[i][j] = y[j]; 
        }
    }

    std::vector<std::vector<double>> p(X, std::vector<double>(Y, 0)), b(X, std::vector<double>(Y, 0));
    b[X / 4][X / 4] = 100;
    b[3 * Y / 4][3 * Y / 4] = -100;

    for (int n = 0; n < T; n++) {
        // Poisson equation calculation
        for (int i = 1; i < X - 1; i++) {
            for (int j = 1; j < Y - 1; j++) {
                p[i][j] = (dy * dy * (p[i][j + 1] + p[i][j - 1]) + dx * dx * (p[i + 1][j] + p[i - 1][j]) - b[i][j] * dy * dy * dx * dx) / (2 * (dy * dy + dx * dx));
            }
        }

        // Apply boundary conditions
        for (int i = 1; i < X - 1; i++) {
            p[i][0] = p[i][1];
            p[i][Y - 1] = p[i][Y - 2];
        }

       const long animatedFig = plt::figure(1);
        plt::ion(); 
        plt::clf();
        plt::plot_surface(nX, nY, p, {}, animatedFig);
        plt::pause(0.05); 
    }

    plt::show(); 
    return 0;
}