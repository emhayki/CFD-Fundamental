
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 9: Laplace equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;
 T = 500;
 X = 31;
 Y = 31;
dx = 2 / (X - 1);
dy = 2 / (Y - 1);

% Create spatial grids
x = linspace(0, 2, X);
y = linspace(0, 2, Y);

% Initial conditions
p = zeros(X, Y);

% Boundary conditions
p(X, :) = y; 

[nX, nY] = meshgrid(x, y);
 
% Main loop for time steps
for n = 1:T
    for i = 2:X-1
        for j = 2:Y-1
           p(i,j) = (dy^2 * (p(i+1,j) + p(i-1,j)) + dx^2 * (p(i,j+1) + p(i,j-1)))/(2 * (dy^2 + dx^2));
        end
    end
    
    % Apply boundary conditions
    p(2:X-1, 1) = p(2:X-1, 2); 
    p(2:X-1, Y) = p(2:X-1, Y-1); 

    % Plot the updated surface
    hFig = figure(1);
    surf(nX, nY, p);
end
