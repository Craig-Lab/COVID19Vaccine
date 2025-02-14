function [fit_parameters,residual,p] = fitting_GCB(param_guess, p, lb, ub)
%setting the optimisation routine specifics
options = optimoptions(@lsqnonlin,'Algorithm', 'trust-region-reflective','MaxFunEval',1000,'DiffMinChange',0.01, 'TolFun', 1e-8, ...  % Function tolerance (smaller is more precise)
    'TolX', 1e-8, ...    % Step size tolerance (smaller is more precise)
   'display','iter-detailed');
% Invoking optimiser
[fit_parameters,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@residualsfunction, param_guess, lb, ub, options);
%------------------------------------------------------------------------
function val = residualsfunction(param)   
    %p.rho_g = param(1);
    p.del_ig= param(1);
    %p.rho_g= param(2);
    p.d_g= param(2);
    tspan = [0 60];
    time_deval = linspace(tspan(1),tspan(2),1e3);
    [time1,sol,solstruc ]= GCBcells(p,tspan);
     yvalstrue= deval(solstruc,time_deval,5);
     [~, idx] = max(yvalstrue);  % Index of the peak value
    peak_time = time_deval(idx);  % Time of the peak
    
    % Objective: Minimize the difference between the peak time and day 8
    val = abs(peak_time - 7);
end
    % val= yvalstrue(idx1)*1e3- dataset; 
end
