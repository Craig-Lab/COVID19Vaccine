%% Objective function to minimize peak time difference
function val = fitting_Bcells(params)
    % Unpack parameters
    p.a = params(1);
    p.del_bt = params(2);
    
    % Time span for ODE simulation
    tspan = [0 60];
    time_deval = linspace(tspan(1), tspan(2), 1e3);
    
    % Run ODE solver
    [~, ~, solstruc] = ActivatedBcells(p, tspan);
    
    % Extract the 4th variable (or the variable of interest)
    yvalstrue = deval(solstruc, time_deval, 4);
    
    % Find the peak value and corresponding time
    [~, idx] = max(yvalstrue);
    peak_time = time_deval(idx);  % Time of the peak
    
    % Objective: Minimize the absolute difference between peak_time and 10 days
    val = abs(peak_time - 10);  % We minimize the absolute difference
end