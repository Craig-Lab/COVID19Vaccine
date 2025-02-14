function [fit_parameters,residual,p] = fitting_T(param_guess, p, lb, ub)
%setting the optimisation routine specifics
options = optimoptions(@lsqnonlin, ...
   'Algorithm', 'levenberg-marquardt', ...
    'MaxFunEvals', 5000, ...               % Increased number of function evaluations
    'MaxIterations', 3000, ...             % Increased number of iterations
    'TolX', 1e-9, ...                      % Tighter tolerance for parameter changes
    'TolFun', 1e-9, ...                    % Tighter tolerance for function values
    'DiffMinChange', 0.001, ...             % Smaller minimum finite difference step
    'display', 'iter-detailed');           % Keep detailed iteration display
% Invoking optimiser
[fit_parameters,~,residual,~,~,~,~] = lsqnonlin(@residualsfunction, param_guess, lb, ub, options);
%------------------------------------------------------------------------
function val = residualsfunction(params)   
    p.del_tv = params(1);
    p.d_t = params(2);
    tspan = [0 60];
    % Run the ODE solver with the updated parameters
    [time, sol, ~] = Tcells(p, tspan);

    % Extract T cell population (Th) from the solution
    Th = sol(3, :);

    % Find the time of the peak for Th
    [~, idx_peak] = max(Th);
    actual_peak_time = time(idx_peak);

    % Target peak time is 6 (assuming 4 days)
    target_peak_time =4;

    % Calculate the squared error (for fmincon optimization)
    val = (actual_peak_time - target_peak_time)^2;

end
    % val= yvalstrue(idx1)*1e3- dataset; 
end
