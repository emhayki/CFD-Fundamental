
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 7: 2D Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

 X = 31;
 Y = 31;
 T = 40;
nu = .05;
dx = 2 / (X - 1);
dy = 2 / (Y - 1);
dt = .25 * dx * dy / nu;

x = linspace(0, 2, X);
y = linspace(0, 2, Y);

u = ones(X,Y);
u(x >= 0.5 & x <= 1, y >= 0.5 & y <= 1) = 2;
 
[nX, nY] = meshgrid(x, y);

for n = 1:T
    un = u;
    for i = 2:X-1
        for j = 2:Y-1
        u(i,j) =  un(i,j) + nu * dt/dx^2 * (un(i+1,j) - 2*un(i,j) + un(i-1,j)) + nu * dt/dy^2 * (un(i,j+1) - 2*un(i,j) + un(i,j-1));
        end
    end

    u(1:Y, 1) = 1; u(1, 1:X) = 1;
    u(1:X, Y) = 1; u(Y, 1:X) = 1;


    % Plot the updated surface
    hFig = figure(1);
    surf(nX, nY, u);

    xlim([0, 2]);
    ylim([0, 2]);
    zlim([1, 2]);
    daspect([1, 1, 1]);
end

