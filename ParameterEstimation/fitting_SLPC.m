function [fit_parameters,residual,p] = fitting_SLPC(param_guess, p, lb, ub)
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
        % Update parameter in structure p
        p.rho_s = param(1);
        p.d_s = param(2);
        % p.beta_s = param(2);  % Uncomment this if you're using beta_s

        % Simulate using SLPC model
        tspan = [0 60];
        [time, sol, ~] = SLPCs(p, tspan);

        % Extract plasmablasts data
        plasmablasts = sol(5, :);  % Assuming sol(5,:) is the plasmablast data
        [~, idx] = max(plasmablasts);  % Find the index of the peak value
        peak_time = time(idx);  % Time of the peak

        % Objective: Minimize the difference between peak time and day 7
        % val should be a vector for lsqnonlin, so we return the scalar difference
        val = abs(peak_time - 8)^2;  % Return the difference from day 7
    end

    % Return the residual as the residual of the fit (in addition to fit parameters)
end