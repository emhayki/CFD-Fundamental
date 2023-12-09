
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6: 2D Convection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

 X = 101;
 Y = 101;
 T = 80;
dx = 2 / (X - 1);
dy = 2 / (Y - 1);
dt = 0.2 * dx;

% Create spatial grids
x = linspace(0, 2, X);
y = linspace(0, 2, Y);

% Preallocate memory
 u = ones(X, Y);
 v = ones(X, Y);
un = ones(X, Y);
vn = ones(X, Y);

u(x >= 0.5 & x <= 1, y >= 0.5 & y <= 1) = 2;
v(x >= 0.5 & x <= 1, y >= 0.5 & y <= 1) = 2;

[nX, nY] = meshgrid(x, y);

% Main loop for time steps
for n = 1:T

    un = u;
    vn = v;

    % Update interior points using the 2D linear convection equation
    for i = 2:X-1
        for j = 2:Y-1
            u(i, j) = un(i, j) - (un(i, j) * dt / dx * (un(i, j) - un(i-1, j))) - (vn(i, j) * (un(i, j) - un(i, j-1)) * dt / dy);
            v(i, j) = vn(i, j) - (un(i, j) * dt / dx * (vn(i, j) - vn(i-1, j))) - (vn(i, j) * (vn(i, j) - vn(i, j-1)) * dt / dy);
        end
    end
    
    % Apply boundary conditions
    u(1:Y, 1) = 1; u(1, 1:X) = 1;
    u(1:X, Y) = 1; u(Y, 1:X) = 1;

    v(1:Y, 1) = 1; v(1, 1:X) = 1;
    v(1:X, Y) = 1; v(Y, 1:X) = 1;

    % Plot the updated surface
    hFig = figure(1);
    surf(nX, nY, u);
    view( -37.5, 30);

end
