
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: 1D Linear Convection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

 X = 41;                    % Number of spatial points
 T = 40;                    % Number of time steps
dt = 0.025;                 % Time step size
dx = 2 / (X - 1);           % Spatial step size
 c = 1;                     % Wave speed

% Initialise spatial grid and initial condition
x = linspace(0, 5, X);
u = ones(1, X);
u(x >= 0.5 & x <= 1) = 2;


for n = 1:T
    un = u;
    for i = 2:X
        u(i) = un(i) - c * (un(i) - un(i-1)) * dt / dx;
    end

    % Plot
    hFig = figure(1);
    set(hFig, 'Position', [360 278 560 249])
    plot(x, u, 'LineWidth', 1);
    axis([0, 2, 0.5, 2.5]);
    title(sprintf('Time Step: %d', n),'fontsize',14);
    set(0,'defaultTextInterpreter','latex');
    xlabel('$x$','fontsize',14);
    ylabel('$u$','fontsize',14);
end

