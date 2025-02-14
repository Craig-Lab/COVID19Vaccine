function [fit_parameters,residual,p] = fitting_dt(dataset, param_guess, p, lb, ub)
%setting the optimisation routine specifics
options = optimoptions(@lsqnonlin, ...
   'Algorithm', 'levenberg-marquardt', ...
    'MaxFunEvals', 5000, ...               % Increased number of function evaluations
    'MaxIterations', 3000, ...             % Increased number of iterations
    'TolX', 1e-9, ...                      % Tighter tolerance for parameter changes
    'TolFun', 1e-9, ...                    % Tighter tolerance for function values
    'DiffMinChange', 0.01, ...             % Smaller minimum finite difference step
    'display', 'iter-detailed');           % Keep detailed iteration display
% Invoking optimiser
[fit_parameters,~,residual,~,~,~,~] = lsqnonlin(@residualsfunction, param_guess, lb, ub, options);
%------------------------------------------------------------------------
function val = residualsfunction(param)          
    p.d_t= param;    % decay rate of T cells
    tspan = [0 60];
     [time,sol,~ ]= Humoral_response_model1(p,tspan);
     yvalstrue= sol(9,:);
     %sol_deval = interp1(time1,sol',time_deval);
     idx = find(time>=30,1);
     val= abs(yvalstrue(idx)- dataset)^2; 
    % val= solstruc- dataset; 
end
%------------------------------------------------------------------------

end