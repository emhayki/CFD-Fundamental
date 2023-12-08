
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: 1D Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close;

 X = 41;                   % Number of spatial points
 T = 20;                   % Number of time steps
 v = 0.1;                  % Diffusion coefficient
dx = 2 / (X - 1);          % Spatial step size
dt = 0.2 * dx^2 / v;       % Time step size


% Initialise spatial grid and initial condition
u = ones(1, X);     
u(1, round(0.5 / dx):round(1 / dx + 1)) = 2;  
x = linspace(0, 2, X);

for n = 1:T
    un = u;
    for i = 2:X-1
        u(i) = v * dt/dx^2 * (un(i+1) - 2*un(i) + un(i-1)) + un(i);
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

