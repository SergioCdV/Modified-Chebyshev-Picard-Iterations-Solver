%% MCPI solver %%
% Author: Sergio Cuevas del Valle
% Date: 19/10/2025

%% Test 1 %% 
% This script discusses the integration of an IVP for the Kepler problem % 

%% Initial data 
N = 10;             % Order of the approximation 
tol = 1E-14;        % Tolerance for convergence

%% Create the solver
mySolver = MCPI( N, tol );

%% Integration 
mu = 1;                     % Gravitational parameter
delta_t = 100;              % Time stan
x0 = [1 0 0 0 1 0];         % Initial conditions
x0 = repmat( x0, N, 1 );    % Initial guess

[t, x, ~] = mySolver.Solve( @(t,s)dynamics(mu, delta_t, t, s), , x0);

%% Results
plot3( x(:,1), x(:,2), x(:,3) )
grid on;
xlabel('x')
xlabel('y')
xlabel('z')

%% Auxiliary functions 
function [ds] = dynamics(mu, delta_t, t, s)
    % Kepler vector field
    ds = delta_t * [s(3:6,:); -mu * s(1:3,:) ./ vecnorm( s(1:3,:) )];
end