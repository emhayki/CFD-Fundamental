#include "matplotlibcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace plt = matplotlibcpp;

////////////////////////////////////////////////////////////
// Step 12: Channel Flow
////////////////////////////////////////////////////////////

int main() {
    const int X = 21;
    const int Y = 21;
    const int T = 50;

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

    const double   F = 1.00;
    const double  dt = 0.01;
    const double  nu = 0.10;
    const double rho = 1.00;

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
    
        // Source term b calculation
        for (int i = 1; i < X - 1; i++) {
            for (int j = 1; j < Y - 1; j++) {
                b[i][j] = rho * (1.0 / dt * ((u[i+1][j] - u[i-1][j]) / (2.0 * dx) + (v[i][j+1] - v[i][j-1]) / (2.0 * dy)) - ((u[i+1][j] - u[i-1][j]) / (2.0 * dx)) * ((u[i+1][j] - u[i-1][j]) / (2.0 * dx)) - 2.0 * ((u[i][j+1] - u[i][j-1]) / (2.0 * dy) * (v[i+1][j] - v[i-1][j]) / (2.0 * dx)) - ((v[i][j+1] - v[i][j-1]) / (2.0 * dy)) * ((v[i][j+1] - v[i][j-1]) / (2.0 * dy)));
            }
        }

        // Boundary Condition
        for (int j = 1; j < Y - 1; j++) b[0][j] = rho * (1.0 / dt * ((u[1][j] - u[X-1][j]) / (2.0 * dx) + (v[0][j+1] - v[0][j-1]) / (2.0 * dy)) - ((u[1][j] - u[X-1][j]) / (2.0 * dx)) * ((u[1][j] - u[X-1][j]) / (2.0 * dx)) - 2.0 * ((u[0][j+1] - u[0][j-1]) / (2.0 * dy) * (v[1][j] - v[X-1][j]) / (2.0 * dx)) - ((v[0][j+1] - v[0][j-1]) / (2.0 * dy)) * ((v[0][j+1] - v[0][j-1]) / (2.0 * dy)));
        for (int j = 1; j < Y - 1; j++) b[X-1][j] = rho * (1.0 / dt * ((u[0][j] - u[X-2][j]) / (2.0 * dx) + (v[X-1][j+1] - v[X-1][j-1]) / (2.0 * dy)) - ((u[0][j] - u[X-2][j]) / (2.0 * dx)) * ((u[0][j] - u[X-2][j]) / (2.0 * dx)) - 2.0 * ((u[X-1][j+1] - u[X-1][j-1]) / (2.0 * dy) * (v[0][j] - v[X-2][j]) / (2.0 * dx)) - ((v[X-1][j+1] - v[X-1][j-1]) / (2.0 * dy)) * ((v[X-1][j+1] - v[X-1][j-1]) / (2.0 * dy)));
        

        // Pressure Poisson equation
        for (int it = 1; it <= 50; it++) {
            std::vector<std::vector<double>> pn = p;
            
            for (int i = 1; i < X - 1; i++) {
                for (int j = 1; j < Y - 1; j++) {
                    p[i][j] = ((pn[i+1][j] + pn[i-1][j]) * dy * dy + (pn[i][j+1] + pn[i][j-1]) * dx * dx - b[i][j] * dx * dx * dy * dy) / (2. * (dx * dx + dy * dy));
                }
            }

            // Boundary Condition
            for (int j = 1; j < Y - 1; j++) p[0][j] = ((pn[1][j] + pn[X-1][j]) * dy * dy + (pn[0][j+1] + pn[0][j-1]) * dx * dx - b[0][j] * dx * dx * dy * dy) / (2. * (dx * dx + dy * dy));
            for (int j = 1; j < Y - 1; j++) p[X-1][j] = ((pn[0][j] + pn[X-2][j]) * dy * dy + (pn[X-1][j+1] + pn[X-1][j-1]) * dx * dx - b[X-1][j] * dx * dx * dy * dy) / (2. * (dx * dx + dy * dy));


            // Wall boundary conditions
            for (int i = 0; i < Y; i++){
                p[i][0] = p[i][1];
                p[i][X-1] = p[i][X-2];
            } 
        }


        std::vector<std::vector<double>> un = u;
        std::vector<std::vector<double>> vn = v;


        // Velocity update
        for (int i = 1; i < X - 1; i++) {
            for (int j = 1; j < Y - 1; j++) {
                u[i][j] = un[i][j] - un[i][j] * dt / dx * (un[i][j] - un[i-1][j]) - vn[i][j] * dt / dy * (un[i][j] - un[i][j-1]) - dt / (2 * rho * dx) * (p[i+1][j] - p[i-1][j]) + nu * dt * ((un[i+1][j] - 2 * un[i][j] + un[i-1][j]) / (dx * dx) + (un[i][j+1] - 2 * un[i][j] + un[i][j-1]) / (dy * dy)) + F * dt;
                v[i][j] = vn[i][j] - un[i][j] * dt / dx * (vn[i][j] - vn[i-1][j]) - vn[i][j] * dt / dy * (vn[i][j] - vn[i][j-1]) - dt / (2 * rho * dy) * (p[i][j+1] - p[i][j-1]) + nu * dt * ((vn[i+1][j] - 2 * vn[i][j] + vn[i-1][j]) / (dx * dx) + (vn[i][j+1] - 2 * vn[i][j] + vn[i][j-1]) / (dy * dy));
            }
        }

        // Boundary Condition
        for (int j = 1; j < Y - 1; j++) {
            u[0][j] = un[0][j] - un[0][j] * dt / dx * (un[0][j] - un[X-1][j]) - vn[0][j] * dt / dy * (un[0][j] - un[0][j-1]) - dt / (2 * rho * dx) * (p[1][j] - p[X-1][j]) + nu * dt * ((un[1][j] - 2 * un[0][j] + un[X-1][j]) / (dx * dx) + (un[0][j+1] - 2 * un[0][j] + un[0][j-1]) / (dy * dy)) + F * dt;
            v[0][j] = vn[0][j] - un[0][j] * dt / dx * (vn[0][j] - vn[X-1][j]) - vn[0][j] * dt / dy * (vn[0][j] - vn[0][j-1]) - dt / (2 * rho * dy) * (p[0][j+1] - p[0][j-1]) + nu * dt * ((vn[1][j] - 2 * vn[0][j] + vn[X-1][j]) / (dx * dx) + (vn[0][j+1] - 2 * vn[0][j] + vn[0][j-1]) / (dy * dy));
        }

        for (int j = 1; j < Y - 1; j++) {
            u[X-1][j] = un[X-1][j] - un[X-1][j] * dt / dx * (un[X-1][j] - un[X-2][j]) - vn[X-1][j] * dt / dy * (un[X-1][j] - un[X-1][j-1]) - dt / (2 * rho * dx) * (p[0][j] - p[X-2][j]) + nu * dt * ((un[0][j] - 2 * un[X-1][j] + un[X-2][j]) / (dx * dx) + (un[X-1][j+1] - 2 * un[X-1][j] + un[X-1][j-1]) / (dy * dy)) + F * dt;
            v[X-1][j] = vn[X-1][j] - un[X-1][j] * dt / dx * (vn[X-1][j] - vn[X-2][j]) - vn[X-1][j] * dt / dy * (vn[X-1][j] - vn[X-1][j-1]) - dt / (2 * rho * dy) * (p[X-1][j+1] - p[X-1][j-1]) + nu * dt * ((vn[0][j] - 2 * vn[X-1][j] + vn[X-2][j]) / (dx * dx) + (vn[X-1][j+1] - 2 * vn[X-1][j] + vn[X-1][j-1]) / (dy * dy));
        }

        for (int i = 0; i < X; i++) {
            u[i][0] = 0; u[i][Y-1] = 0;
            v[i][0] = 0; v[i][Y-1] = 0;
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