% PARAMETER BASELINE VALUES
load_parameters_new

params =p;

% Parameter Labels 
    PRCC_var={'\delta_{LV}','d_l','d_V','\delta_{TV}','d_T','\delta_{BT}','\rho_G','\rho_S', 'd_B','\delta_{IG}','\beta_{G}',...
    'p_G','d_G','p_P','d_P','d_S','d_M','p_M','\beta_M','p_P2','ht','\alpha_P','\alpha_S','d_A', '\rho_I','d_I','SI'};% checked!

%% TIME SPAN OF THE SIMULATION
t_end=60; % length of the simulations
tspan=linspace(0, 60,1e3);   % time points where the output is calculated
time_points=[round(t_end/2) t_end]; % time points of interest for the US analysis

% INITIAL CONDITION FOR THE ODE MODEL       
%y0 = [0.16; 1e-10; 1e-10; 0; 0; 0.015; 1e-7; 5e-6];
%params.IC = [1;0;0;0;0;0;0;0;0;0];

% Variables Labels
y_var_label={'L';'V';'Th';'B';'GB';'LP';'SP';'M';'A';'I'};