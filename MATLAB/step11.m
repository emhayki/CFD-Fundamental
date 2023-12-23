
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 11: Cavity Flow 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

X = 41;
Y = 41;
T = 100;

  c = 1.000;
 dt = 0.001;
 nu = 0.100;
rho = 1.000;

dx = 2 / (X - 1);
dy = 2 / (Y - 1);

x = linspace(0, 2, X);
y = linspace(0, 2, Y);

[nX, nY] = meshgrid(x, y);

b = zeros(X, Y);
p = zeros(X, Y);
u = zeros(X, Y);
v = zeros(X, Y);

% Lid-driven cavity condition
u(end, :) = 1;

for n = 1:T

    un = u;
    vn = v;

    % Source term b calculation
    for i = 2:X-1
        for j = 2:Y-1
            b(i,j) = rho * (1 / dt * ((u(i+1,j) - u(i-1,j)) / (2 * dx) + (v(i,j+1) - v(i,j-1)) / (2 * dy)) - ((u(i+1,j) - u(i-1,j)) / (2 * dx))^2 - 2 * ((u(i,j+1) - u(i,j-1)) / (2 * dy) * (v(i+1,j) - v(i-1,j)) / (2 * dx)) - ((v(i,j+1) - v(i,j-1)) / (2 * dy))^2);
        end
    end

    % Pressure Poisson equation
    for it = 1:50
        pn = p;
        for i = 2:X-1
            for j = 2:Y-1
                p(i,j) = ((pn(i+1,j) + pn(i-1,j)) * dy^2 + (pn(i,j+1) + pn(i,j-1)) * dx^2 - b(i,j) * dx^2 * dy^2) / (2 * (dx^2 + dy^2));
            end
        end

        % Pressure boundary conditions
        p(:,1) = p(:,2);  
        p(1,:) = p(2,:);   
        p(end,:) = 0;     
        p(:,end) = p(:,end-1); 
    end

  
    % Velocity update
    for i = 2:X-1
        for j = 2:Y-1
            u(i,j) = un(i,j) - un(i,j) * dt / dx * (un(i,j) - un(i-1,j)) - vn(i,j) * dt / dy * (un(i,j) - un(i,j-1)) - dt / (2 * rho * dx) * (p(i+1,j) - p(i-1,j)) + nu * dt * ((un(i+1,j) - 2 * un(i,j) + un(i-1,j)) / dx^2 + (un(i,j+1) - 2 * un(i,j) + un(i,j-1)) / dy^2);
            v(i,j) = vn(i,j) - un(i,j) * dt / dx * (vn(i,j) - vn(i-1,j)) - vn(i,j) * dt / dy * (vn(i,j) - vn(i,j-1)) - dt / (2 * rho * dy) * (p(i,j+1) - p(i,j-1)) + nu * dt * ((vn(i+1,j) - 2 * vn(i,j) + vn(i-1,j)) / dx^2 + (vn(i,j+1) - 2 * vn(i,j) + vn(i,j-1)) / dy^2);
        end
    end

    % Velocity boundary conditions
    u(1,:) = 0; u(:,1) = 0; u(:,end) = 1; u(end,:) = 0;
    v(1,:) = 0; v(:,1) = 0; v(:,end) = 0; v(end,:) = 0;

  
    hFig = figure(1);
    contourf(nX, nY, p);
    hold on;
    quiver(nX,nY,u',v',3);
    colorbar
    xlabel('X');
    ylabel('Y');
    title('Pressure and Velocity Field');

end

