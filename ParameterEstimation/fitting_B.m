function [fit_parameters,residual,p] = fitting_B(param_guess, p, lb, ub)
%setting the optimisation routine specifics
options = optimoptions(@lsqnonlin, ...
    'Algorithm', 'Levenberg-Marquardt', ...
    'MaxFunEvals', 5000, ...               % Increased number of function evaluations
    'MaxIterations', 3000, ...             % Increased number of iterations
    'TolX', 1e-9, ...                      % Tighter tolerance for parameter changes
    'TolFun', 1e-9, ...                    % Tighter tolerance for function values
    'DiffMinChange', 0.001, ...             % Smaller minimum finite difference step
    'display', 'iter-detailed');           % Keep detailed iteration display
% Invoking optimiser
[fit_parameters,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@residualsfunction, param_guess, lb, ub, options);
%------------------------------------------------------------------------
function val = residualsfunction(param)   
    p.a = param(1);
    p.del_bt = param(1);
    tspan = [0 60];
    time_deval = linspace(tspan(1),tspan(2),1e3);
    [time1,sol,solstruc ]= ActivatedBcells(p,tspan);
     yvalstrue= deval(solstruc,time_deval,4);
     [ymax, idx] = max(yvalstrue);  % Index of the peak value
    peak_time = time_deval(idx);  % Time of the peak
    
    % Objective: Minimize the difference between the peak time and day 10
    val = abs(peak_time - 5)^2;

end
    % val= yvalstrue(idx1)*1e3- dataset; 
end
