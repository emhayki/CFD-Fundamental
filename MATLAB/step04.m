
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: 1D Burgers' Equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

X = 101;                      % Number of spatial points
T = 100;                      % Number of time steps
v = 0.07;                     % Diffusion coefficient
dx = 2 * pi / (X - 1);        % Spatial step size
dt = dx * v;                  % Time step size

% Initialize spatial grid and initial condition
 x = linspace(0, 2*pi, X);
 u = zeros(1, X); 
uA = zeros(1, X); 

% Initialise initial conditions using the analytical solution
for i = 1:X
     u(i) = AnalyticalSolution(0, v, x(i));
    uA(i) = AnalyticalSolution(100*dt, v, x(i));
end



for n = 1:T
    un = u;
    
    % Discretisation of the Burgers' equation
    for i = 2:X-1
        u(i) = (v * dt/dx^2 * (un(i+1) - 2*un(i) + un(i-1))) - (un(i) * dt/dx * (un(i) - un(i-1))) + un(i);
    end
    
    % Boundary conditions
    u(1) = (v * dt/dx^2 * (un(2) - 2*un(1) + un(X-1))) - (un(1) * dt/dx * (un(1) - un(X-1))) + un(1);
    u(X) = u(1);

    % Plot the solution at each time step
    hFig = figure(1);
    set(hFig, 'Position', [360 278 560 249])
    plot(x, u, 'LineWidth', 1);
    axis([0, 6, 0.5, 8]);
    title(sprintf('Time Step: %d', n),'fontsize',14);
    set(0,'defaultTextInterpreter','latex');
    xlabel('$x$','fontsize',14);
    ylabel('$u$','fontsize',14);
end

% Plot the analytical solution
hold on
plot(x, uA, 'LineWidth', 1);
legend("Computational", "Analytical")

% Analytical solution function
function u = AnalyticalSolution(t, v, x)
     phi = exp(-(4*t - x)^2/(4*v*(1 + t))) + exp(-(2*pi + 4*t - x)^2/(4*v*(1 + t)));
    dphi = (exp(-(2*pi + 4*t - x)^2/(4*v*(1 + t)))*(4*pi + 8*t - 2*x))/(4*v*(1 + t)) + (exp(-(4*t - x)^2/(4*v*(1 + t)))*(8*t - 2*x))/(4*v*(1 + t));
       u = 4 - (2*dphi*v)/phi;
end
