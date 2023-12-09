
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: 2D Linear Convection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

 X = 81;
 Y = 81;
 T = 100;
 c = 1;
dx = 2 / (X - 1);
dy = 2 / (Y - 1);
dt = 0.2 * dx;

% Create spatial grids
x = linspace(0, 2, X);
y = linspace(0, 2, Y);
u = ones(X, Y);
u(x >= 0.5 & x <= 1, y >= 0.5 & y <= 1) = 2;

[nX, nY] = meshgrid(x, y);

% Main loop for time steps
for n = 1:T
    un = u;

    % Update interior points using the 2D linear convection equation
    for i = 2:X-1
        for j = 2:Y-1
            u(i, j) = un(i, j) - c * (un(i, j) - un(i-1, j)) * dt / dx - c * (un(i, j) - un(i, j-1)) * dt / dx;
        end
    end
    
    % Apply boundary conditions
    u(1:Y, 1) = 1; u(1, 1:X) = 1;
    u(1:X, Y) = 1; u(Y, 1:X) = 1;

    % Plot the updated surface
    hFig = figure(1);
    surf(nX, nY, u, 'FaceColor', 'interp');
    title(sprintf('Time Step: %d', n), 'fontsize', 14);
    grid on;
    colormap(jet); 
end
