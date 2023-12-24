
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 12: Channel Flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

% Parameters
X = 21;
Y = 21;
T = 50;

% Domain and grid spacing
xmin = 0;    xmax = 2;
ymin = 0;    ymax = 1;

dx = (xmax - xmin) / (X - 1);
dy = (ymax - ymin) / (Y - 1);

% Create spatial grids
x = linspace(xmin, xmax, X);
y = linspace(ymin, ymax, Y);

% Flow parameters
  F = 1.00;
 dt = 0.01;
 nu = 0.10;
rho = 1.00;

% Meshgrid for plotting
[nX, nY] = meshgrid(x, y);

b = zeros(X, Y);
p = zeros(X, Y);
u = zeros(X, Y);
v = zeros(X, Y);

pn = zeros(X, Y);
un = zeros(X, Y);
vn = zeros(X, Y);

for n = 1:T

    % Source term b calculation
    for i = 2:X-1
        for j = 2:Y-1
            b(i, j) = rho * (1 / dt * ((u(i+1, j) - u(i-1, j)) / (2 * dx) + (v(i, j+1) - v(i, j-1)) / (2 * dy)) - ((u(i+1, j) - u(i-1, j)) / (2 * dx))^2 - 2 * ((u(i, j+1) - u(i, j-1)) / (2 * dy) * (v(i+1, j) - v(i-1, j)) / (2 * dx)) - ((v(i, j+1) - v(i, j-1)) / (2 * dy))^2);
        end
    end

    % Periodic Boundary Condition @ x = 0
    for j = 2:Y-1
        b(1, j) = rho * (1 / dt * ((u(2, j) - u(X, j)) / (2 * dx) + (v(1, j+1) - v(1, j-1)) / (2 * dy)) - ((u(2, j) - u(X, j)) / (2 * dx))^2 - 2 * ((u(1, j+1) - u(1, j-1)) / (2 * dy) * (v(2, j) - v(X, j)) / (2 * dx)) - ((v(1, j+1) - v(1, j-1)) / (2 * dy))^2);
    end

    % Periodic Boundary Condition @ x = 2
    for j = 2:Y-1
        b(X, j) = rho * (1 / dt * ((u(1, j) - u(X-1, j)) / (2 * dx) + (v(X, j+1) - v(X, j-1)) / (2 * dy)) - ((u(1, j) - u(X-1, j)) / (2 * dx))^2 - 2 * ((u(X, j+1) - u(X, j-1)) / (2 * dy) * (v(1, j) - v(X-1, j)) / (2 * dx)) - ((v(X, j+1) - v(X, j-1)) / (2 * dy))^2);
    end

    % Pressure Poisson equation
    for it = 1:50
        pn = p;
        for i = 2:X-1
            for j = 2:Y-1
                p(i, j) = ((pn(i+1, j) + pn(i-1, j)) * dy^2 + (pn(i, j+1) + pn(i, j-1)) * dx^2 - b(i, j) * dx^2 * dy^2) / (2 * (dx^2 + dy^2));
            end
        end

        % Periodic Boundary Condition @ x = 0
        for j = 2:Y-1
            p(1, j) = ((pn(2, j) + pn(X, j)) * dy^2 + (pn(1, j+1) + pn(1, j-1)) * dx^2 - b(1, j) * dx^2 * dy^2) / (2 * (dx^2 + dy^2));
        end

        % Periodic Boundary Condition @ x = 2
        for j = 2:Y-1
            p(X, j) = ((pn(1, j) + pn(X-1, j)) * dy^2 + (pn(X, j+1) + pn(X, j-1)) * dx^2 - b(X, j) * dx^2 * dy^2) / (2 * (dx^2 + dy^2));
        end

        % Wall boundary conditions
        p(:, 1) = p(:, 2);
        p(:, end) = p(:, end-1);
    end

    un = u;
    vn = v;

    % Velocity update
    for i = 2:X-1
        for j = 2:Y-1
            u(i, j) = un(i, j) - un(i, j) * dt / dx * (un(i, j) - un(i-1, j)) - vn(i, j) * dt / dy * (un(i, j) - un(i, j-1)) - dt / (2 * rho * dx) * (p(i+1, j) - p(i-1, j)) + nu * dt * ((un(i+1, j) - 2 * un(i, j) + un(i-1, j)) / dx^2 + (un(i, j+1) - 2 * un(i, j) + un(i, j-1)) / dy^2) + F * dt;
            v(i, j) = vn(i, j) - un(i, j) * dt / dx * (vn(i, j) - vn(i-1, j)) - vn(i, j) * dt / dy * (vn(i, j) - vn(i, j-1)) - dt / (2 * rho * dy) * (p(i, j+1) - p(i, j-1)) + nu * dt * ((vn(i+1, j) - 2 * vn(i, j) + vn(i-1, j)) / dx^2 + (vn(i, j+1) - 2 * vn(i, j) + vn(i, j-1)) / dy^2);
        end
    end

    % Periodic Boundary Condition @ x = 0
    for j = 2:Y-1
        u(1, j) = un(1, j) - un(1, j) * dt / dx * (un(1, j) - un(X, j)) - vn(1, j) * dt / dy * (un(1, j) - un(1, j-1)) - dt / (2 * rho * dx) * (p(2, j) - p(X, j)) + nu * dt * ((un(2, j) - 2 * un(1, j) + un(X, j)) / dx^2 + (un(1, j+1) - 2 * un(1, j) + un(1, j-1)) / dy^2) + F * dt;
        v(1, j) = vn(1, j) - un(1, j) * dt / dx * (vn(1, j) - vn(X, j)) - vn(1, j) * dt / dy * (vn(1, j) - vn(1, j-1)) - dt / (2 * rho * dy) * (p(1, j+1) - p(1, j-1)) + nu * dt * ((vn(2, j) - 2 * vn(1, j) + vn(X, j)) / dx^2 + (vn(1, j+1) - 2 * vn(1, j) + vn(1, j-1)) / dy^2);
    end

    % Periodic Boundary Condition @ x = 2
    for j = 2:Y-1
        u(X, j) = un(X, j) - un(X, j) * dt / dx * (un(X, j) - un(X-1, j)) - vn(X, j) * dt / dy * (un(X, j) - un(X, j-1)) - dt / (2 * rho * dx) * (p(1, j) - p(X-1, j)) + nu * dt * ((un(1, j) - 2 * un(X, j) + un(X-1, j)) / dx^2 + (un(X, j+1) - 2 * un(X, j) + un(X, j-1)) / dy^2) + F * dt;
        v(X, j) = vn(X, j) - un(X, j) * dt / dx * (vn(X, j) - vn(X-1, j)) - vn(X, j) * dt / dy * (vn(X, j) - vn(X, j-1)) - dt / (2 * rho * dy) * (p(X, j+1) - p(X, j-1)) + nu * dt * ((vn(1, j) - 2 * vn(X, j) + vn(X-1, j)) / dx^2 + (vn(X, j+1) - 2 * vn(X, j) + vn(X, j-1)) / dy^2);
    end

    % Velocity boundary conditions
    u(:, 1) = 0; u(:, Y) = 0;
    v(:, 1) = 0; v(:, Y) = 0;

    % Plot velocity vectors
    hFig = figure(2);
    quiver(nX, nY, u', v', 1);
end
