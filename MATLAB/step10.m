
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 10: 2D Poisson Equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

 X = 50;
 Y = 50;
 T = 100;

xmin = 0;    xmax = 2;
ymin = 0;    ymax = 1;

dx = (xmax - xmin)/(X - 1);
dy = (ymax - ymin)/(Y - 1);

% Create spatial grids
x = linspace(xmin, xmax, X);
y = linspace(xmin, ymax, Y);

% Initial conditions
p = zeros(X, Y);
b = zeros(X, Y);

b(floor(1/4 * X), floor(1/4 * X)) =  100;
b(floor(3/4 * Y), floor(3/4 * Y)) = -100;


[nX, nY] = meshgrid(x, y);

% Main loop for time steps
for n = 1:T
    for i = 2:X-1
        for j = 2:Y-1
           p(i,j) = (dy^2 * (p(i,j+1) + p(i,j-1)) + dx^2 * (p(i+1,j) + p(i-1,j)) - b(i,j)*dy^2*dx^2)/(2 * (dy^2 + dx^2));
        end
    end
    
    % Apply boundary conditions
    p(2:X-1, 1) = p(2:X-1, 2); 
    p(2:X-1, Y) = p(2:X-1, Y-1); 

    % Plot the updated surface
    hFig = figure(1);
    surf(nX, nY, p);
    zlim([-0.05 0.05])
end
