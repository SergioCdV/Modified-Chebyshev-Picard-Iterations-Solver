%% Modified Chebyshev Picard Iteration Solver %% 
% Sergio Cuevas del Valle
% Date: 19/10/25
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Solve %%
% This function implements the iteration procedure to solve a given BP

function [x, Output] = Solve(obj, dynamics, problemType, x0)
    % Pre-allocation 
    m = size(x0,2);                         % Dimension of the state space
    chi = zeros( obj.N, m );                % Boundary conditions term
    X = zeros(m * obj.N, obj.maxIter);      % Evolution of the trajectory

    % Set up of the method 
    GoOn = true;                            % Convergence boolean
    iter = 1;                               % Initial iteration
    error = 1E3;                            % Initial error

    switch problemType
        case BP.IVP
            hndl_ = @(x) ( [2 * x(1,:); zeros(1, size(x,2))] );

        case BP.TBVP
            hndl_ = @(x) ( [x(1,:) + x(end,:); 0.5 * ( x(end,:) - x(1,:) ) ] );
    end

    % Main loop
    while ( GoOn && iter <= obj.maxIter )
        % Initial iteration 
        X(:,iter) = reshape( x0, [], 1 );

        % Evaluate the dynamics 
        g = dynamics(tau, x0); 

        % Compute the beta coefficients as a function of the problem
        chi(1:2,:) = hndl_( x );
        beta = obj.Ca * g + chi; 

        % Compute the new state trajectory
        x = obj.Cx * beta;

        % Convergence check 
        dX = x - x0; 
        dX = sqrt( dot(dX, dX, 2) );

        if ( error < tol && norm(dX) < tol )
            GoOn = false;
        else
            % Update the iteration scheme
            x0 = x;

            error = norm(dX);
            iter = iter + 1;
        end
    end

    % Final output 
    Output.Result     = ~GoOn; 
    Output.Iterations = iter; 
    Output.Error      = error; 
    Output.Evolution  = X(:,1:iter);
end