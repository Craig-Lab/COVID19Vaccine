%% HCW
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
load('hcw_ant_ind.mat')

%%
Nboosters = 6;%number of boosters
time_initial = time;
hcw_initial = H;
hcw_booster1 = zeros(10,4000+1000*Nboosters,length(hcw_dt));
time_booster1 = zeros(1,4000+1000*Nboosters);
booster = 1;
lag = 360;

for i = 1:length(hcw_dt)
    tic
        p.lam2 = hcw_lam2(i);
        p.lam3 =hcw_lam3(i);
        p.d_t  = hcw_dt(i);
        hcw_sol = squeeze(hcw_initial(:,:,i));
        hcw_time = time_initial;
    for j = 1:Nboosters %number of boosters 
        p.L0 = hcw_sol(1,end)+booster;
        p.V0 = hcw_sol(2,end);
        p.Th0 = hcw_sol(3,end);
        p.B0 =  hcw_sol(4,end);
        p.GB0 =  hcw_sol(5,end);
        p.LP0 =  hcw_sol(6,end);
        p.SP0 = hcw_sol(7,end);
        p.M0= hcw_sol(8,end);
        p.A0= hcw_sol(9,end);
        p.I0 =  hcw_sol(10,end);
        tspan0 = [hcw_time(end) hcw_time(end)+lag];
        [sol, time,solstruc] = Humoral_response_model3(p, tspan0);
        time_deval = linspace(tspan0(1),tspan0(2),1e3);
        sol_deval = deval(solstruc,time_deval); 
        hcw_sol = [hcw_sol, sol_deval];
        hcw_time = [hcw_time, time_deval];
    end
    hcw_booster1(:,:,i) = hcw_sol;
    time_booster1 =hcw_time;
i
end

save('hcw_booster1_yearly.mat',"time_booster1","hcw_booster1")


%
clear all
format long
hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
load('hcw_ant_ind.mat')
Nboosters = 6;%number of boosters
time_initial = time;
hcw_initial = H;
hcw_booster05 = zeros(10,4000+1000*Nboosters,length(hcw_dt));
time_booster05 = zeros(1,4000+1000*Nboosters);
booster = 0.5;
lag = 360;

for i = 1:length(hcw_dt)
    tic
        p.lam2 = hcw_lam2(i);
        p.lam3 =hcw_lam3(i);
        p.d_t  = hcw_dt(i);
        hcw_sol = squeeze(hcw_initial(:,:,i));
        hcw_time = time_initial;
    for j = 1:Nboosters %number of boosters 
        p.L0 = hcw_sol(1,end)+booster;
        p.V0 = hcw_sol(2,end);
        p.Th0 = hcw_sol(3,end);
        p.B0 =  hcw_sol(4,end);
        p.GB0 =  hcw_sol(5,end);
        p.LP0 =  hcw_sol(6,end);
        p.SP0 = hcw_sol(7,end);
        p.M0= hcw_sol(8,end);
        p.A0= hcw_sol(9,end);
        p.I0 =  hcw_sol(10,end);
        tspan0 = [hcw_time(end) hcw_time(end)+lag];
        [sol, time,solstruc] = Humoral_response_model3(p, tspan0);
        time_deval = linspace(tspan0(1),tspan0(2),1e3);
        sol_deval = deval(solstruc,time_deval); 
        hcw_sol = [hcw_sol, sol_deval];
        hcw_time = [hcw_time, time_deval];
    end
    hcw_booster05(:,:,i) = hcw_sol;
    time_booster05 =hcw_time;
i
end

save('hcw_booster05_yearly.mat',"time_booster05","hcw_booster05")


%% senior
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
load('senior_ant_ind.mat')
Nboosters = 6;%number of boosters
time_initial = time;
senior_initial = H;
senior_booster1 = zeros(10,4000+1000*Nboosters,length(senior_dt));
time_booster1 = zeros(1,4000+1000*Nboosters);
booster = 1;
lag = 360;

for i = 1:length(senior_dt)
    tic
        p.lam2 = senior_lam2(i);
        p.lam3 =senior_lam3(i);
        p.d_t  = senior_dt(i);
        senior_sol = squeeze(senior_initial(:,:,i));
        senior_time = time_initial;
    for j = 1:Nboosters %number of boosters 
        p.L0 = senior_sol(1,end)+booster;
        p.V0 = senior_sol(2,end);
        p.Th0 = senior_sol(3,end);
        p.B0 =  senior_sol(4,end);
        p.GB0 =  senior_sol(5,end);
        p.LP0 =  senior_sol(6,end);
        p.SP0 = senior_sol(7,end);
        p.M0= senior_sol(8,end);
        p.A0= senior_sol(9,end);
        p.I0 =  senior_sol(10,end);
        tspan0 = [senior_time(end) senior_time(end)+lag];
        [sol, time,solstruc] = Humoral_response_model3(p, tspan0);
        time_deval = linspace(tspan0(1),tspan0(2),1e3);
        sol_deval = deval(solstruc,time_deval); 
        senior_sol = [senior_sol, sol_deval];
        senior_time1 = [senior_time, time_deval];
    end
    senior_booster1(:,:,i) = senior_sol;
    time_booster1 =senior_time;

end

save('senior_booster1_yearly.mat',"time_booster1","senior_booster1")
%% senior half dose
clear all
format long
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
load('senior_ant_ind.mat')
Nboosters = 6;%number of boosters
time_initial = time;
senior_initial = H;
senior_booster05 = zeros(10,4000+1000*Nboosters,length(senior_dt));
time_booster05 = zeros(1,4000+1000*Nboosters);
booster = 0.5;
lag = 360;

for i = 1:length(senior_dt)
    tic
        p.lam2 = senior_lam2(i);
        p.lam3 =senior_lam3(i);
        p.d_t  = senior_dt(i);
        senior_sol = squeeze(senior_initial(:,:,i));
        senior_time = time_initial;
    for j = 1:Nboosters %number of boosters 
        p.L0 = senior_sol(1,end)+booster;
        p.V0 = senior_sol(2,end);
        p.Th0 = senior_sol(3,end);
        p.B0 =  senior_sol(4,end);
        p.GB0 =  senior_sol(5,end);
        p.LP0 =  senior_sol(6,end);
        p.SP0 = senior_sol(7,end);
        p.M0= senior_sol(8,end);
        p.A0= senior_sol(9,end);
        p.I0 =  senior_sol(10,end);
        tspan0 = [senior_time(end) senior_time(end)+lag];
        [sol, time,solstruc] = Humoral_response_model3(p, tspan0);
        time_deval = linspace(tspan0(1),tspan0(2),1e3);
        sol_deval = deval(solstruc,time_deval); 
        senior_sol = [senior_sol, sol_deval];
        senior_time = [senior_time, time_deval];
    end

    i
    senior_booster05(:,:,i) = senior_sol;
    time_booster05 =senior_time;

end

save('senior_booster05_yearly.mat',"time_booster05","senior_booster05")
%
%senior half-dose 6m
clear all
format long
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
load('senior_ant_ind.mat')
Nboosters = 12;%number of boosters
time_initial = time;
senior_initial = H;
senior_booster05_6m = zeros(10,4000+1000*Nboosters,length(senior_dt));
time_booster05_6m = zeros(1,4000+1000*Nboosters);
booster = 0.5;
lag = 180;

for i = 1:length(senior_dt)
    tic
        p.lam2 = senior_lam2(i);
        p.lam3 =senior_lam3(i);
        p.d_t  = senior_dt(i);
        senior_sol = squeeze(senior_initial(:,:,i));
        senior_time = time_initial;
    for j = 1:Nboosters %number of boosters 
        p.L0 = senior_sol(1,end)+booster;
        p.V0 = senior_sol(2,end);
        p.Th0 = senior_sol(3,end);
        p.B0 =  senior_sol(4,end);
        p.GB0 =  senior_sol(5,end);
        p.LP0 =  senior_sol(6,end);
        p.SP0 = senior_sol(7,end);
        p.M0= senior_sol(8,end);
        p.A0= senior_sol(9,end);
        p.I0 =  senior_sol(10,end);
        tspan0 = [senior_time(end) senior_time(end)+lag];
        [sol, time,solstruc] = Humoral_response_model3(p, tspan0);
        time_deval = linspace(tspan0(1),tspan0(2),1e3);
        sol_deval = deval(solstruc,time_deval); 
        senior_sol = [senior_sol, sol_deval];
        senior_time = [senior_time, time_deval];
    end
    i
    senior_booster05_6m(:,:,i) = senior_sol;
    time_booster05_6m =senior_time;

end

save('senior_booster05_6m.mat',"time_booster05_6m","senior_booster05_6m")

%%
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
load('hcw_ant_ind.mat')
Nboosters = 8;%number of boosters
time_initial = time;
hcw_initial = H;
hcw_booster05_6m = zeros(10,4000+1000*Nboosters,length(hcw_dt));
time_booster05_6m = zeros(1,4000+1000*Nboosters);
booster = 0.5;
lag = 180;

for i = 1:length(hcw_dt)
    tic
        p.lam2 = hcw_lam2(i);
        p.lam3 =hcw_lam3(i);
        p.d_t  = hcw_dt(i);
        hcw_sol = squeeze(hcw_initial(:,:,i));
        hcw_time = time_initial;
    for j = 1:Nboosters %number of boosters 
        p.L0 = hcw_sol(1,end)+booster;
        p.V0 = hcw_sol(2,end);
        p.Th0 = hcw_sol(3,end);
        p.B0 =  hcw_sol(4,end);
        p.GB0 =  hcw_sol(5,end);
        p.LP0 =  hcw_sol(6,end);
        p.SP0 = hcw_sol(7,end);
        p.M0= hcw_sol(8,end);
        p.A0= hcw_sol(9,end);
        p.I0 =  hcw_sol(10,end);
        tspan0 = [hcw_time(end) hcw_time(end)+lag];
        [sol, time,solstruc] = Humoral_response_model3(p, tspan0);
        time_deval = linspace(tspan0(1),tspan0(2),1e3);
        sol_deval = deval(solstruc,time_deval); 
        hcw_sol = [hcw_sol, sol_deval];
        hcw_time = [hcw_time, time_deval];
    end
    hcw_booster05_6m(:,:,i) = hcw_sol;
    time_booster05_6m =hcw_time;
end
save('hcw_booster05_6m.mat',"time_booster05_6m","hcw_booster05_6m")
