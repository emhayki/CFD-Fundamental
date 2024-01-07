#include "matplotlibcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace plt = matplotlibcpp;

////////////////////////////////////////////////////////////
// Step 9: Laplace equation
////////////////////////////////////////////////////////////

int main() {
    // Define simulation parameters
    const int X = 31;                           // Number of points along X-axis
    const int Y = 31;                           // Number of points along Y-axis
    const int T = 500;                          // Total number of time steps

    const double dx = 2. / (X - 1);             // Step size in the X direction
    const double dy = 2. / (Y - 1);             // Step size in the Y direction

    // Create spatial grids
    std::vector<double> x(X), y(Y);

    for (int i = 0; i < X; i++)
        x[i] = (2 * i) / (X - 1.0);

    for (int i = 0; i < Y; i++)
        y[i] = (2 * i) / (Y - 1.0);

    std::vector<std::vector<double>> nX(X, std::vector<double>(Y)), nY(X, std::vector<double>(Y));
    for (int i = 0; i < X; ++i) {
        for (int j = 0; j < Y; ++j) {
            nX[i][j] = x[i]; 
            nY[i][j] = y[j]; 
        }
    }

    std::vector<std::vector<double>> p(X, std::vector<double>(Y));

    for (int i = 0; i < X; i++) p[X-1][i] = y[i];
 
    // Time-stepping loop
    for (int n = 0; n < T; n++) {
         for (int i = 1; i < X - 1; ++i) {
            for (int j = 1; j < Y - 1; ++j) {
                p[i][j] = (dy * dy * (p[i+1][j]  + p[i-1][j] ) + dx * dx * (p[i][j+1]  + p[i][j-1] ))/(2 * (dy * dy  + dx * dx));
            }
        }

        // Apply boundary conditions
        for (int i = 1; i < X - 1; ++i) {
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