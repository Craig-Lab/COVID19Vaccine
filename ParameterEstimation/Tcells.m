function [time, sol, solstruc] = Tcells(p, tspan)

    % Initial conditions [Lipid nanoparticles (L), Vaccinated cells (V), Helper T cells (Th)]
    % Use values from 'p' if they exist; otherwise, fall back on the default [1;0;0]
    if isfield(p, 'L0') && isfield(p, 'V0') && isfield(p, 'Th0')
        p.IC = [p.L0; p.V0; p.Th0];
    else
        p.IC = [1; 0; 0];  % Default values if not provided
    end

    % ODE solver options
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 1e-2);

    % Solve the system of ODEs using ode45
    solstruc = ode45(@(t, y) odefun(t, y, p), tspan, p.IC, opts);

    % Extract solution (time points and corresponding solution values)
    time = solstruc.x;        % Time points where the solution was computed
    sol = solstruc.y;         % Solution values for [L; V; Th] at those time points

    % If needed, you can interpolate the solution at specific time points
    % interpolated_time = linspace(tspan(1), tspan(end), 1000);
    % sol = deval(solstruc, interpolated_time);

    %------------------------------------------------------------------------
    % ODE function describing the system
    function dydt = odefun(t, y, p)
        L = y(1);   % Lipid nanoparticles
        V = y(2);   % Vaccinated cells
        Th = y(3);  % Helper T cells

        % Define the rate equations using parameters in 'p'
        
        dL = -p.del_lv * L - p.d_l * L;
        dV = p.del_lv * L - p.d_v * V;
        dTh = p.del_tv * V - p.d_t * Th;

        % Return the derivatives as a column vector
        dydt = [dL; dV; dTh];
    end

end