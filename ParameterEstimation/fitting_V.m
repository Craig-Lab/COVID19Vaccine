function [fit_parameters,residual,p] = fitting_V(param_guess, p, lb, ub)
%setting the optimisation routine specifics
options = optimoptions(@lsqnonlin, ...
   'Algorithm', 'levenberg-marquardt', ...
    'MaxFunEvals', 5000, ...               % Increased number of function evaluations
    'MaxIterations', 3000, ...             % Increased number of iterations
    'TolX', 1e-9, ...                      % Tighter tolerance for parameter changes
    'TolFun', 1e-9, ...                    % Tighter tolerance for function values
    'DiffMinChange', 0.1, ...             % Smaller minimum finite difference step
    'display', 'iter-detailed');           % Keep detailed iteration display
% Invoking optimiser
[fit_parameters,~,residual,~,~,~,~] = lsqnonlin(@residualsfunction, param_guess, lb, ub, options);
%------------------------------------------------------------------------
function val = residualsfunction(param)   
    p.d_v = param(1);
   % p.d_l = param(2);

    tspan = [0 60];
    [time1,sol,~ ]= Vcells(p,tspan);
    vaccinatedcells = sol(2,:);
    idx = find(time1 > 28, 1);

    [~, idx_peak] = max(vaccinatedcells);
    actual_peak_time = time1(idx_peak);

    target_peak_time = 2;

    % Calculate the squared error (for fmincon optimization)

    val = abs(vaccinatedcells(idx) - 1e-6)^2+ (actual_peak_time - target_peak_time)^2;
end
    % val= yvalstrue(idx1)*1e3- dataset; 
end
