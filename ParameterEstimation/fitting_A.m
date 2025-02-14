function [fit_parameters,residual,p] = fitting_A(param_guess, p, lb, ub)
%setting the optimisation routine specifics
options = optimoptions(@lsqnonlin, ...
    'Algorithm', 'levenberg-marquardt', ...      % Try this algorithm for potentially faster convergence
    'MaxFunctionEvaluations', 2000, ...
    'MaxIterations', 1000, ...
    'TolX', 1e-6, ...
    'DiffMinChange',1e-6, ...
    'FiniteDifferenceStepSize', 0.1, ...
    'Display', 'iter-detailed');
% Invoking optimiser
[fit_parameters,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@residualsfunction, param_guess, lb, ub, options);
%------------------------------------------------------------------------pl
function val = residualsfunction(param)   
    p.alpha_p = param(1);
    p.alpha_s= param(2);

    tspan = [0 60];
   
    [time,sol,~]= Antibodies(p,tspan);
     antibodies = sol(9,:);
    
     idx = find(time>=30,1); % Index of the peak value
     
    anti = antibodies(idx);  % Time of the peak
    [~, idx_peak] = max(antibodies);
    actual_peak_time = time(idx_peak);

    % Target peak time is 6 (assuming 6 days)
    target_peak_time = 21;
%abs(anti - 88.86)^2+
    % Objective: Minimize the difference between the peak time and day 10
   %val = abs(actual_peak_time-target_peak_time)^2;
    val = abs(anti - 88.86)^2+abs(actual_peak_time-target_peak_time)^2;
end
    % val= yvalstrue(idx1)*1e3- dataset; 
end
