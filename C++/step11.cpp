#include "matplotlibcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace plt = matplotlibcpp;

////////////////////////////////////////////////////////////
// Step 11: Cavity Flow 
////////////////////////////////////////////////////////////

int main() {
    int X = 41;
    int Y = 41;
    int T = 100;

    double dt = 0.001;
    double nu = 0.100;
    double rho = 1.000;

    double dx = 2.0 / (X - 1);
    double dy = 2.0 / (Y - 1);

    std::vector<std::vector<double>> nX(X, std::vector<double>(Y));
    std::vector<std::vector<double>> nY(X, std::vector<double>(Y));
    std::vector<std::vector<double>>  b(X, std::vector<double>(Y));
    std::vector<std::vector<double>>  p(X, std::vector<double>(Y));
    std::vector<std::vector<double>>  u(X, std::vector<double>(Y));
    std::vector<std::vector<double>>  v(X, std::vector<double>(Y));

    // Initialize nX and nY
    for (int i = 0; i < X; i++) {
        for (int j = 0; j < Y; j++) {
            nX[i][j] = 2.0 * i / (X - 1);
            nY[i][j] = 2.0 * j / (Y - 1);
        }
    }

    // Lid-driven cavity condition
    for (int i = 0; i < Y; i++) {
        u[X - 1][i] = 1.0;
    }


    for (int n = 0; n < T; n++) {
        std::vector<std::vector<double>> un = u;
        std::vector<std::vector<double>> vn = v;

        // Source term b calculation
        for (int i = 1; i < X - 1; i++) {
            for (int j = 1; j < Y - 1; j++) {
                b[i][j] = rho * (1.0 / dt * ((u[i+1][j] - u[i-1][j]) / (2. * dx) + (v[i][j+1] - v[i][j-1]) / (2. * dy)) - ((u[i+1][j] - u[i-1][j]) / (2. * dx)) * ((u[i+1][j] - u[i-1][j]) / (2. * dx)) - 2. * ((u[i][j+1] - u[i][j-1]) / (2. * dy) * (v[i+1][j] - v[i-1][j]) / (2. * dx)) - ((v[i][j+1] - v[i][j-1]) / (2. * dy)) * ((v[i][j+1] - v[i][j-1]) / (2. * dy)));
            }
        }

        // Pressure Poisson equation
        for (int it = 1; it <= 50; it++) {
            std::vector<std::vector<double>> pn = p;
            
            for (int i = 1; i < X - 1; i++) {
                for (int j = 1; j < Y - 1; j++) {
                    p[i][j] = ((pn[i+1][j] + pn[i-1][j]) * dy * dy + (pn[i][j+1] + pn[i][j-1]) * dx * dx - b[i][j] * dx * dx * dy * dy) / (2. * (dx * dx + dy * dy));
                }
            }

            // Pressure boundary conditions
            for (int i = 0; i < X; i++) {
                p[i][0] = p[i][1];
                p[i][Y-1] = p[i][Y-2];
            }

            for (int j = 0; j < Y; j++) {
                p[0][j] = p[1][j];
                p[X-1][j] = 0.0;
            }
        }

        // Velocity update
        for (int i = 1; i < X - 1; i++) {
            for (int j = 1; j < Y - 1; j++) {
                u[i][j] = un[i][j] - un[i][j] * dt / dx * (un[i][j] - un[i-1][j]) - vn[i][j] * dt / dy * (un[i][j] - un[i][j-1]) - dt / (2. * rho * dx) * (p[i+1][j] - p[i-1][j]) + nu * dt * ((un[i+1][j] - 2. * un[i][j] + un[i-1][j]) / (dx * dx) + (un[i][j+1] - 2. * un[i][j] + un[i][j-1]) / (dy * dy));
                v[i][j] = vn[i][j] - un[i][j] * dt / dx * (vn[i][j] - vn[i-1][j]) - vn[i][j] * dt / dy * (vn[i][j] - vn[i][j-1]) - dt / (2. * rho * dy) * (p[i][j+1] - p[i][j-1]) + nu * dt * ((vn[i+1][j] - 2. * vn[i][j] + vn[i-1][j]) / (dx * dx) + (vn[i][j+1] - 2. * vn[i][j] + vn[i][j-1]) / (dy * dy));
            }
        }

        // Velocity boundary conditions
        for (int i = 0; i < X; i++) {
            u[i][0] = 0.0;
            v[i][0] = 0.0;
            u[i][Y-1] = 1.0;
            v[i][Y-1] = 0.0;
        }

        for (int i = 0; i < Y; i++) {
            u[0][i] = 0.0;
            v[0][i] = 0.0;
            u[X-1][i] = 0.0;
            v[X-1][i] = 0.0;
        }


    plt::figure(1); 
    std::vector<double> X1, Y1, P1, U1, V1;

    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            X1.push_back(j * dx);
            Y1.push_back(i * dy);
            P1.push_back(p[i][j]);
            U1.push_back(u[j][i]);
            V1.push_back(v[j][i]);
        }
    }

    const long animatedFig = plt::figure(1);
    plt::ion(); 
    plt::clf();
    plt::show();
    plt::quiver(X1, Y1, U1, V1);
    plt::pause(0.1);
    }


    return 0;
}