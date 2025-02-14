function [fit_parameters,residual,p] = fitting_BM(param_guess, p, lb, ub)
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
%------------------------------------------------------------------------
function val = residualsfunction(param)   
    p.p_m= param(1);
    tspan2 = [0 180];
    p.lam1 = 0;
    p.lam2 = 0;
    [time,sol,~] = Humoral_response_model2(p,tspan2);
    GCB = sol(5,:);
    
     idx = find(time>=30,1); % Index of the peak value
     
    GCB30 = GCB(idx);  % Time of the peak
   
    val = abs(GCB30  - 1e-6)^2;
end
    % val= yvalstrue(idx1)*1e3- dataset; 
end
