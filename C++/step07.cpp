#include "matplotlibcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace plt = matplotlibcpp;

////////////////////////////////////////////////////////////
// Step 7: 2D Diffusion
////////////////////////////////////////////////////////////

int main() {
    // Define simulation parameters
    const int X = 31;                        // Number of points along X-axis
    const int Y = 31;                        // Number of points along Y-axis
    const int T = 40;                        // Total number of time steps

    const double nu = .05;           
    const double dx = 2. / (X - 1);           // Step size in the X direction
    const double dy = 2. / (Y - 1);           // Step size in the Y direction
    const double dt = .25 * dx * dy / nu;     // Time step size

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

    std::vector<std::vector<double>> u(X, std::vector<double>(Y)), un(X, std::vector<double>(Y));

    for (int i = 0; i < X; i++){
        for (int j = 0; j < Y; j++)
            u[i][j] = ((x[i] >= 0.5 && x[i] <= 1) && (y[j] >= 0.5 && y[j] <= 1)) ? 2.0 : 1.0;     
    }
 
    // Time-stepping loop
    for (int n = 0; n < T; n++) {
        un = u;
        
        for (int i = 1; i < X-1; i++) {
            for (int j = 1; j < Y-1; j++)
            u[i][j] =  un[i][j] + nu * dt/(dx*dx) * (un[i+1][j] - 2*un[i][j] + un[i-1][j]) + nu * dt/(dy*dy) * (un[i][j+1] - 2*un[i][j]+ un[i][j-1]);
        }

        // Boundary conditions
        for (int i = 0; i < X; i++) u[i][0] = 1.;
        for (int i = 0; i < Y; i++) u[0][i] = 1.;

        for (int i = 0; i < X; i++) u[i][Y-1] = 1.;
        for (int i = 0; i < Y; i++) u[X-1][i] = 1.;
 

    const long animatedFig = plt::figure(1);
    plt::ion(); 
    plt::clf();
    plt::plot_surface(nX, nY, u, {}, animatedFig);
    plt::pause(0.1); 
    }
    
    plt::show(); 
    return 0;
}
