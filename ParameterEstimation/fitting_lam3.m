function [fit_parameters,residual,p] = fitting_lam3( dataset, param_guess, p, lb, ub)
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
    lam3= param;    % decay rate of T cells
   % p.IC = [p.L0;p.V0;p.Th0;p.B0;p.P0;p.M0;p.A0;p.I0];
   tspan = [0 360];
    % call the solver, and convert to log10
     [time,sol,~]= Humoral_response_model3(p,tspan,p.lam2,lam3);
     yvalstrue= sol(9,:);
     idx1 = find(time>=30,1);
     idx2 = find(time>=90,1);
     idx3 = find(time>=180,1);
     idxx = [idx1,idx2,idx3];
     val= abs(yvalstrue(idxx)*1e3- dataset).^2; 
end
%------------------------------------------------------------------------

end