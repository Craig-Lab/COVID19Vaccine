function [fit_parameters,residual,p] = fitting_LLPCs(param_guess, p, lb, ub)
%setting the optimisation routine specifics
options = optimoptions(@lsqnonlin,'Algorithm', 'trust-region-reflective','MaxFunEval',1000,'DiffMinChange',0.1, 'TolFun', 1e-8, ...  % Function tolerance (smaller is more precise)
    'TolX', 1e-8, ...    % Step size tolerance (smaller is more precise)
   'display','iter-detailed');
% Invoking optimiser% Invoking optimiser
[fit_parameters,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@residualsfunction, param_guess, lb, ub, options);
%------------------------------------------------------------------------
function val = residualsfunction(param)   
    %p.rho_g = param(1);
    p.p_p= param(1);
    tspan = [0 60];
    [time,sol,~ ]= LLPCs(p,tspan);
    idx = find(time>=10,1);
    llpc= sol(7,1:idx);
    s= sol(6,1:idx);
    tt= time(1:idx);
    total_plasmablasts = trapz(tt,s);
    total_llpcs = trapz(tt, llpc);
   % [~, idx] = max(s);  % Index of the peak value
 
    % Objective: Minimize the difference between the peak time and day 10
    val = abs(1.5*total_llpcs- total_plasmablasts)^2;
end
    % val= yvalstrue(idx1)*1e3- dataset; 
end