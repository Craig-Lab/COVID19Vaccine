%% fit d_V
clear all
rng(1)
p.d_l = 0.1;
p.del_lv =0.15;% 0.15; %0.246 -0.1 active translation time 4 days, I suppose (d_L+del_LV) depends on this time.

% Number of iterations
num_trials = 1;
prevs_trials = 0; % The number of existing trials of the same type

for current_trial = 1:num_trials
    a_lb = 0.1;
    a_ub =1;
    a_guess = (a_ub-a_lb).*rand(1,1) + a_lb;
    par_guess = a_guess;

    lb = a_lb;
    ub = a_ub;
    trial = prevs_trials+current_trial;

    % Run the optimizer
    [fit_parameters, residual, p] = fitting_V(par_guess, p, lb, ub); 
    p.d_v = fit_parameters(1);
   % p.d_l = fit_parameters(2);
end
    tspan1 = [0 60];
    [time1, sol1,solstruc1]= Vcells(p,tspan1);
    time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
    sol1_deval= deval(solstruc1,time_deval1);
   
    figure(1)
    subplot(2,1,1)
    plot(time1,sol1(1,:))
    subplot(2,1,2)
    plot(time1,sol1(2,:))
%% fit parameters of T follicular helper T cells:
clear all
rng(123)
p.d_l = 0.1;
p.del_lv = 0.15; %0.246 -0.1
p.d_v = 0.89;%updated


num_trials = 1;
prevs_trials = 0; % The number of existing trials of the same type

for current_trial = 1:num_trials
    a_lb = 0.1;
    a_ub =1;
    a_guess = (a_ub-a_lb).*rand(1,1) + a_lb;

    del_lb = 0.4;
    del_ub = 1;
    del_guess = (del_ub-del_lb).*rand(1,1) + del_lb;
    
    
    par_guess = [a_guess,del_guess];
    lb = [a_lb,del_lb];
    ub = [a_ub,del_ub];

    trial = prevs_trials+current_trial;

    [fit_parameters, residual, p] = fitting_T(par_guess, p, lb, ub); 
   % p.d_v = fit_parameters(1);
    p.del_tv = fit_parameters(1);
    p.d_t= fit_parameters(2);
end

    tspan1 = [0 60];
    [time1, sol1,solstruc1]= Tcells(p,tspan1);
    time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
    sol1_deval= deval(solstruc1,time_deval1);
    [a1, b1]=max(sol1_deval(2,:));
    [a2, b2]=max(sol1_deval(3,:));
    time_deval1(b1)
     time_deval1(b2)
    %
    figure(1)
    subplot(3,1,1)
    plot(time1,sol1(1,:))
    subplot(3,1,2)
    plot(time1,sol1(2,:)*1e6)
    subplot(3,1,3)   
    plot(time1,sol1(3,:)*1e6)
%% fit parameters of activated B cells:
clear all
rng(123)
p.d_l = 0.1;
p.del_lv = 0.15; %0.246 -0.1
p.d_v = 0.89;%updated
p.del_tv = 0.73;
p.d_t= 0.57;
%p.d_b = 0.02; %0.115;%half-life 6 days
p.d_b = 0.02;
% Number of iterations

num_trials = 1;
prevs_trials = 0; % The number of existing trials of the same type

for current_trial = 1:num_trials
    a_lb = 0.8;
    a_ub = 1.5;
    a_guess = (a_ub-a_lb).*rand(1,1) + a_lb;

    del_lb = 0.1;
    del_ub = 1;
    del_guess = (del_ub-del_lb).*rand(1,1) + del_lb;

    par_guess = [a_guess, del_guess];
    lb = [a_lb, del_lb];
    ub = [a_ub, del_ub];

    trial = prevs_trials+current_trial;
    % Run the optimizer
    [fit_parameters, residual, p] = fitting_B(par_guess, p, lb, ub); 
    p.a = fit_parameters(1);
    p.del_bt = fit_parameters(2);
end

% Final fitted parameters
disp(['  a = ', num2str(p.a)]);
disp(['  del_bt = ', num2str(p.del_bt)]);

    p.a = fit_parameters(1);
    p.del_bt= fit_parameters(2);
    tspan1 = [0 60];
    [time1, sol1,solstruc1]= ActivatedBcells(p,tspan1);
    time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
    sol1_deval= deval(solstruc1,time_deval1);
    [a1, b1]=max(sol1_deval(3,:));
    [a2, b2]=max(sol1_deval(4,:));
    time_deval1(b1)
    time_deval1(b2)
    figure(1)
    subplot(4,1,1)
    plot(time1,sol1(1,:))
    subplot(4,1,2)
    plot(time1,sol1(2,:)*1e6)
    subplot(4,1,3)   
    plot(time1,sol1(3,:)*1e6)
    subplot(4,1,4)   
    plot(time1,sol1(4,:)*1e6)
 %% fitting SLPBs
clear all
rng(1)
p.d_l = 0.1;
p.del_lv = 0.15; %0.246 -0.1
p.d_v = 0.89;%updated
p.del_tv = 0.73;
p.d_t= 0.57;
p.d_b = 0.02;
p.sum = 1.29;
p.del_bt =0.36;

num_trials = 1;
prevs_trials = 0; % The number of existing trials of the same type

for current_trial = 1:num_trials
    a_lb = 0.1;
    a_ub = p.sum;
   
    % Generate random initial guesses within bounds
    a_guess = a_lb + (a_ub - a_lb) * rand;

    del_lb = 0.2;
    del_ub = 0.33;
    del_guess = (del_ub-del_lb).*rand(1,1) + del_lb;

    par_guess = [a_guess, del_guess];
    lb = [a_lb, del_lb];
    ub = [a_ub, del_ub];


    trial = prevs_trials+current_trial;

    [fit_parameters, residual, p] = fitting_SLPC(par_guess, p, lb, ub); 
    p.rho_s = fit_parameters(1);
    p.d_s = fit_parameters(2);


end


tspan1 = [0 60];
[time1, sol1,solstruc1]= SLPCs(p,tspan1);
time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
sol1_deval= deval(solstruc1,time_deval1);
[~, b2]=max(sol1_deval(5,:));
time_deval1(b2)

figure(1)
subplot(5,1,1);
plot(time1,sol1(1,:))
%xlim([0 5])
subplot(5,1,2);
plot(time1,sol1(2,:))
%xlim([0 5])
subplot(5,1,3);
plot(time1,sol1(3,:))
%xlim([0 5])
title('Th')
subplot(5,1,4);
plot(time1,sol1(4,:))
title('B')
subplot(5,1,5);
plot(time1,sol1(5,:))
title('SLPCs')
%% fit GC Bcells
clear all
rng(1)
p.d_l = 0.1;
p.del_lv = 0.15; %0.246 -0.1
p.d_v = 0.89;%updated
p.del_tv = 0.73;
p.d_t= 0.57;
p.d_b = 0.02; %0.115;%half-life 6 days
p.sum = 1.29;
p.del_bt =0.36;
p.d_s = 0.29;%lifespan 3-5 days(0.33, 0.2) 
p.rho_s =0.6;


p.SI= 0.03;%4.5;%
p.ht = 1;%%%%%
p.rho_i =2.07;
p.d_i = 3.47;%half-life 5hours

p.beta_g = 2.77;
p.p_g =0.1;
p.rho_g = p.sum-p.rho_s;

tspan = [0 60];  

% Number of iterations
num_trials = 1;
prevs_trials = 0; % The number of existing trials of the same type

for current_trial = 1:num_trials
    a_lb = 0.1;
    a_ub = 1;
   
    del_lb = 0.5;
    del_ub = 0.8;
    % Generate random initial guesses within bounds
 
    a_guess = (a_ub - a_lb) .* rand(1,1) + a_lb;
    del_guess = (del_ub - del_lb) .* rand(1,1) + del_lb;
    par_guess = [a_guess, del_guess];
    lb = [a_lb, del_lb];
   ub = [a_ub, del_ub];
    %par_guess = a_guess;
   % lb = a_lb;
   %ub = a_ub;
    trial = prevs_trials+current_trial;

    [fit_parameters, residual, p] = fitting_GCB(par_guess, p, lb, ub); 
   % p.rho_g = fit_parameters(1);
    p.del_ig =fit_parameters(1);
   % p.rho_g = fit_parameters(2);
    p.d_g = fit_parameters(2);

end

    tspan1 = [0 60];
    [time1, sol1,solstruc1]= GCBcells(p,tspan1);
    time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
    sol1_deval= deval(solstruc1,time_deval1);
    [a1, b1]=max(sol1_deval(4,:));
    [a2, b2]=max(sol1_deval(5,:));
    time_deval1(b1)
    time_deval1(b2)
    %
figure(1)
subplot(7,1,1);
plot(time1,sol1(1,:))
%xlim([0 5])
subplot(7,1,2);
plot(time1,sol1(2,:))
%xlim([0 5])
subplot(7,1,3);
plot(time1,sol1(3,:))
%xlim([0 5])
title('Th')
subplot(7,1,4);
plot(time1,sol1(4,:))
%xlim([0 5])
title('B')
subplot(7,1,5);
plot(time1,sol1(5,:))
title('GCBs')
subplot(7,1,6);
plot(time1,sol1(6,:))
title('SLPCs')
subplot(7,1,7);
plot(time1,sol1(7,:))
title('IL-21')
%% fit LLPCs
clear all
rng(1)
p.d_l = 0.1;
p.del_lv = 0.15; %0.246 -0.1
p.d_v = 0.89;%updated
p.del_tv = 0.73;
p.d_t= 0.57;
p.d_b = 0.02; %half-life of 5 weeks
p.sum = 1.29;
p.del_bt =0.36;
p.d_s = 0.29;%lifespan 3-5 days(0.33, 0.2) 
p.rho_s =0.6;


p.SI= 0.03;%4.5;%
p.ht = 1;%%%%%
p.rho_i =2.07;
p.d_i = 3.47;%half-life 5hours

p.beta_g = 2.77;
p.p_g =0.1;
p.rho_g = p.sum-p.rho_s;



tspan = [0 60];  

%p.del_ig =1.37;
p.del_ig =0.72;
p.d_g =0.48;
p.d_p =  0.0112;% 


% Number of iterations
num_trials = 1;
prevs_trials = 0; % The number of existing trials of the same type
for current_trial = 1:num_trials
    a_lb = 0.1;
    a_ub = 0.5;

    % Generate random initial guesses within bounds
 
    a_guess = (a_ub - a_lb) .* rand(1,1) + a_lb;
    par_guess = a_guess;
    lb = a_lb;
    ub = a_ub;
    trial = prevs_trials+current_trial;

    [fit_parameters, residual, p] = fitting_LLPCs(par_guess, p, lb, ub); 
   % p.rho_g = fit_parameters(1);
    p.p_p=fit_parameters(1);

end

tspan1 = [0 60];
[time1, sol1,solstruc1]= LLPCs(p,tspan1);
time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
sol1_deval= deval(solstruc1,time_deval1);
[a1, b1]=max(sol1_deval(4,:));
[a2, b2]=max(sol1_deval(5,:));
time_deval1(b1)
time_deval1(b2)
figure(1)
subplot(6,1,1);
plot(time1,sol1(1,:))
%xlim([0 5])
subplot(6,1,2);
plot(time1,sol1(2,:)*1e9)
%xlim([0 5])
subplot(6,1,3);
plot(time1,sol1(3,:)*1e9)
%xlim([0 5])
title('Th')
subplot(6,1,4);
plot(time1,sol1(4,:)*1e9)
%xlim([0 5])
title('B')
subplot(6,1,5);
plot(time1,sol1(5,:)*1e9)
title('GCBs')
subplot(6,1,6);
plot(time1,sol1(6,:)*1e9)
hold on
plot(time1,sol1(7,:)*1e9)
title('LLPCs')
%% fit antibodies
rng(1)
p.d_l = 0.1;
p.del_lv = 0.15; %0.246 -0.1
p.d_v = 0.89;%updated
p.del_tv = 0.73;
p.d_t= 0.57;
p.d_b = 0.02; %0.115;%half-life 6 days
p.sum = 1.29;
p.del_bt =0.36;
p.d_s = 0.29;%lifespan 3-5 days(0.33, 0.2) 
p.rho_s =0.6;


p.SI= 0.03;%4.5;%
p.ht = 1;%%%%%
p.rho_i =2.07;
p.d_i = 3.47;%half-life 5hours

p.beta_g = 2.77;
p.p_g =0.1;
p.rho_g = p.sum-p.rho_s;



tspan = [0 60];  

%p.del_ig =1.37;
p.del_ig =0.72;
p.d_g =0.48;
p.d_p =  0.0112;% 

p.p_p = 0.26;
p.d_m = 0.0029;
p.beta_m = 0.009;%fitted 0.022;%0.022;%22;%0.022;%0.022*100;% 0.0266;
p.d_a =  0.033; %halflife 60 days; mean 90
    
% Number of iterations
num_trials = 1;
prevs_trials = 0; % The number of existing trials of the same type

for current_trial = 1:num_trials
    a_lb = 10;
    a_ub = 30;
   
    del_lb = 200;
    del_ub = 300;
    % Generate random initial guesses within bounds
 
    a_guess = (a_ub - a_lb) *rand+ a_lb;
    del_guess = (del_ub - del_lb) *rand + del_lb;
    par_guess = [a_guess, del_guess];
    lb = [a_lb, del_lb];
    ub = [a_ub, del_ub];
    %par_guess = a_guess;
    %lb = a_lb;
    %ub = a_ub;
    trial = prevs_trials+current_trial;

    [fit_parameters, residual, p] = fitting_A(par_guess, p, lb, ub); 
   % p.rho_g = fit_parameters(1);
    p.alpha_p =fit_parameters(1);
    %p.p_g = fit_parameters(2);
    p.alpha_s = fit_parameters(2);

end
%
    tspan1 = [0 60];
    % [time1, sol1,solstruc1]= LLPCs(p,tspan1);
    [time1, sol1,solstruc1]= Antibodies(p,tspan1);
    time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
    sol1_deval= deval(solstruc1,time_deval1);

figure(1)
subplot(8,1,1);
plot(time1,sol1(1,:))
%xlim([0 5])
subplot(8,1,2);
plot(time1,sol1(2,:))
%xlim([0 5])
subplot(8,1,3);
plot(time1,sol1(3,:))
%xlim([0 5])
title('Th')
subplot(8,1,4);
plot(time1,sol1(4,:))
%xlim([0 5])
title('B')
subplot(8,1,5);
plot(time1,sol1(5,:))
title('GCBs')
subplot(8,1,6);
plot(time1,sol1(6,:))
hold on
plot(time1,sol1(7,:))
title('LLPCs')
subplot(8,1,7);
plot(time1,sol1(8,:))
title('MBCs')
subplot(8,1,8);
plot(time1,sol1(9,:))
hold on
scatter(30, 88.86)
title('Antibodies')

%%
clear all
load_parameters_new
tspan1 = [0 60];
[time1, sol1,solstruc1]= Humoral_response_model1(p,tspan1);
time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
sol1_deval= deval(solstruc1,time_deval1);
figure(1)
subplot(8,1,1);
plot(time1,sol1(1,:)*1e3)
%xlim([0 5])
subplot(8,1,2);
plot(time1,sol1(2,:)*1e9)
%xlim([0 5])
subplot(8,1,3);
plot(time1,sol1(3,:)*1e9)
%xlim([0 5])
title('Th')
subplot(8,1,4);
plot(time1,sol1(4,:)*1e9)
%xlim([0 5])
title('B')
subplot(8,1,5);
plot(time1,sol1(5,:)*1e9)
title('GCBs')
subplot(8,1,6);
plot(time1,sol1(6,:)*1e9)
hold on
plot(time1,sol1(7,:)*1e9)
title('LLPCs')

subplot(8,1,7);
plot(time1,sol1(8,:)*1e9)
title('MBCs')
subplot(8,1,8);
plot(time1,sol1(9,:))
hold on
scatter(30, 88.86)
title('Antibodies')
%%
p.lam1=0;
p.lam2=0;
p.L0 = sol1_deval(1,1000)+1;
p.V0 = sol1_deval(2,1000);
p.Th0 = sol1_deval(3,1000);
p.B0 =  sol1_deval(4,1000);
p.GB0 =  sol1_deval(5,1000);
p.LP0 =  sol1_deval(6,1000);
p.SP0 = sol1_deval(7,1000);
p.M0= sol1_deval(8,1000);
p.A0= sol1_deval(9,1000);
p.I0 =  sol1_deval(10,1000);
p.p_p2 = 0.9;
p.beta_m = 0.009;%% leterature
%p.p_p2 = 0.95;
rng(123)
% Number of iterations
num_trials = 1;
prevs_trials = 0; % The number of existing trials of the same type
for current_trial = 1:num_trials
    a_lb = 0.5;
    a_ub =1;

    % Generate random initial guesses within bounds
 
    a_guess = (a_ub - a_lb) .* rand(1,1) + a_lb;
    par_guess = a_guess;
    lb = a_lb;
    ub = a_ub;
    trial = prevs_trials+current_trial;

    [fit_parameters, residual, p] = fitting_BM(par_guess, p, lb, ub); 
   % p.rho_g = fit_parameters(1);
   p.p_m =fit_parameters(1);

end
tspan2 = [0 180];
%p.IC = [p.L0;p.V0;p.Th0;p.B0;p.P0;p.M0;p.A0;p.I0];
[time2,sol2,solstruc2] = Humoral_response_model2(p,tspan2);
time_deval2 = linspace(tspan2(1),tspan2(2),1e3);
sol2_deval = deval(solstruc2,time_deval2);
GCB = sol2(5,:);

idx = find(time2>=30,1); % Index of the peak value

GCB30 = GCB(idx)  % Time of the peak
%
figure(2)
subplot(8,1,1);
plot(time2,sol2(1,:)*1e3)
%xlim([0 5])
subplot(8,1,2);
plot(time2,sol2(2,:)*1e9)
%xlim([0 5])
subplot(8,1,3);
plot(time2,sol2(3,:)*1e9)
%xlim([0 5])
title('Th')
subplot(8,1,4);
plot(time2,sol2(4,:)*1e9)
%xlim([0 5])
title('B')
subplot(8,1,5);
plot(time2,sol2(5,:)*1e9)
title('GCBs')
subplot(8,1,6);
plot(time2,sol2(6,:)*1e9)
title('LLPCs')
subplot(8,1,7);
plot(time2,sol2(8,:)*1e9)
title('MBCs')
subplot(8,1,8);
plot(time2,sol2(9,:))
title('Antibodies')

%%

clear all
load_parameters_new
tspan1 = [0 60];
[time1, sol1,solstruc1]= Humoral_response_model1(p,tspan1);
time_deval1 = linspace(tspan1(1),tspan1(2),1e3);
sol1_deval= deval(solstruc1,time_deval1);
p.lam2=7.6;
p.lam3=0;
p.L0 = sol1_deval(1,1000)+1;
p.V0 = sol1_deval(2,1000);
p.Th0 = sol1_deval(3,1000);
p.B0 =  sol1_deval(4,1000);
p.GB0 =  sol1_deval(5,1000);
p.LP0 =  sol1_deval(6,1000);
p.SP0 = sol1_deval(7,1000);
p.M0= sol1_deval(8,1000);
p.A0= sol1_deval(9,1000);
p.I0 =  sol1_deval(10,1000);
tspan2 = [0 180];
[time2,sol2,solstruc2] = Humoral_response_model2(p,tspan2);
time_deval2 = linspace(tspan2(1),tspan2(2),1e3);
sol2_deval = deval(solstruc2,time_deval2);
p.lam2=7.6;
p.lam3=5.7;
p.L0 = sol2_deval(1,1000)+1;
p.V0 = sol2_deval(2,1000);
p.Th0 = sol2_deval(3,1000);
p.B0 =  sol2_deval(4,1000);
p.GB0 =  sol2_deval(5,1000);
p.LP0 =  sol2_deval(6,1000);
p.SP0 = sol2_deval(7,1000);
p.M0= sol2_deval(8,1000);
p.A0= sol2_deval(9,1000);
p.I0 =  sol2_deval(10,1000);
tspan3 = [0 360];
[time3,sol3,solstruc3] = Humoral_response_model2(p,tspan3);
time_deval3 = linspace(tspan3(1),tspan3(2),1e3);
sol3_deval = deval(solstruc3,time_deval3);
GCB = sol2(5,:);
time_total = [time_deval1,time_deval2+60, time_deval3+60+180];
sol_total = [sol1_deval,sol2_deval,sol3_deval];
figure(2)
subplot(8,1,1);
plot(time_total,sol_total(1,:)*1e3)
%xlim([0 5])
subplot(8,1,2);
plot(time_total,sol_total(2,:)*1e9)
%xlim([0 5])
subplot(8,1,3);
plot(time_total,sol_total(3,:)*1e9)
%xlim([0 5])
title('Th')
subplot(8,1,4);
plot(time_total,sol_total(4,:)*1e9)
%xlim([0 5])
title('B')
subplot(8,1,5);
plot(time_total,sol_total(5,:)*1e9)
title('GCBs')
subplot(8,1,6);
plot(time_total,sol_total(6,:)*1e9)
title('LLPCs')
subplot(8,1,7);
plot(time_total,sol_total(8,:)*1e9)
title('MBCs')
subplot(8,1,8);
plot(time_total,sol_total(9,:))
title('Antibodies')
%%
figure(3)
subplot(10,1,1);
plot(time1,sol1(1,:))
hold on
plot(time2,sol2(1,:))
%xlim([0 5])
xlim([0 60])
subplot(10,1,2);
plot(time1,sol1(2,:))
hold on
plot(time2,sol2(2,:))
xlim([0 60])
%xlim([0 5])
subplot(10,1,3);
plot(time1,sol1(3,:))
hold on
plot(time2,sol2(3,:))
xlim([0 60])
%xlim([0 5])
title('Th')
subplot(10,1,4);
plot(time1,sol1(4,:))
hold on
plot(time2,sol2(4,:))
xlim([0 60])
%xlim([0 5])
title('B')
subplot(10,1,5);
plot(time1,sol1(5,:))
hold on
plot(time2,sol2(5,:))
xlim([0 60])
title('GCBs')
subplot(10,1,6);
plot(time1,sol1(6,:))
hold on
plot(time2,sol2(6,:))
xlim([0 60])
title('LLPCs')
subplot(10,1,7);
plot(time1,sol1(7,:))
hold on
plot(time2,sol2(7,:))
xlim([0 60])
title('SLPCs')
subplot(10,1,8);
plot(time1,sol1(8,:))
hold on
plot(time2,sol2(8,:))
xlim([0 60])
title('MBCs')
subplot(10,1,9);
plot(time1,sol1(9,:))
hold on
plot(time2,sol2(9,:))
xlim([0 60])
title('Antibodies')
subplot(10,1,10);
plot(time1,sol1(10,:))
hold on
plot(time2,sol2(10,:))
xlim([0 60])
title('IL-21')
%%
t = [time1, 60+time2];
S = [sol1, sol2];
figure(4)
plot(t,p.alpha_p*S(6,:),'r' )
hold on
plot(t,p.alpha_s*S(7,:),'k'  )
%%
figure(4)
plot(time2, (1-p.p_g)*p.p_p*p.beta_g*sol2(5,:))
hold on
plot(time2, p.beta_m*(1-p.p_m)*p.p_p*sol2(8,:))
