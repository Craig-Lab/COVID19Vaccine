clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
dose3 =0.5;


H  = zeros(10,4000,length(hcw_lam2));
for i = 1:length(hcw_lam2)
    p.lam2 = hcw_lam2(i);
    p.lam3=hcw_lam3(i);
    p.d_t  = hcw_dt(i);

tic
 [sol, time] = Model_3doses(p,dose3);
toc
H(:,:,i)= sol;
i
end


save('hcw_ant_ind.mat','H','time')

%% senior
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
dose3 =1;


H  = zeros(10,4000,length(senior_lam2));
for i = 1:length(senior_lam2)
    p.lam2 = senior_lam2(i);
    p.lam3=senior_lam3(i);
    p.d_t  = senior_dt(i);

tic
 [sol, time] = Model_3doses(p,dose3);
toc
H(:,:,i)= sol;
i
end

save('senior_ant_ind.mat','H','time')


%%
hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
% lipid nanoparticles
load_parameters_new()
dose3 =1;


H  = zeros(10,4000,length(hcw_lam2));
for i = 1:length(hcw_lam2)
    p.lam2 = hcw_lam2(i);
    p.lam3=hcw_lam3(i);
    p.d_t  = hcw_dt(i);

tic
 [sol, time] = Model_3doses(p,dose3);
toc
H(:,:,i)= sol;
i
end


save('hcw_ant_ind_fulldose.mat','H','time')

%% Generating neutralization curves for seniors
hcw_par_data = readtable('hcw_neu_par.xlsx');
hcw_par = table2array(hcw_par_data(:,2:4));

load('hcw_ant_ind.mat');
LH =squeeze(H(9,:,:))'./1e3;
HH_dose = zeros(length(hcw_par),4000);
for i = 1:length(hcw_par)
    HH_dose(i,:) =real(hcw_par(i,1)+(1-hcw_par(i,1))* LH(i,:).^hcw_par(i,3)./(LH(i,:).^hcw_par(i,3)+hcw_par(i,2)^hcw_par(i,3)))*100;
end
HH = HH_dose;
Time = time;
save('hcw_neu_ind.mat',"HH",'Time')

%% Generating neutralization curves for seniors
senior_par_data = readtable('senior_neu_par.xlsx');
senior_par = table2array(senior_par_data(:,2:4));
load('senior_ant_ind.mat');
LS =squeeze(H(9,:,:))'./1e3;
SS_dose = zeros(length(senior_par),4000);
for i = 1:length(senior_par)
    SS_dose(i,:) =real(senior_par(i,1)+(1-senior_par(i,1))* LS(i,:).^senior_par(i,3)./(LS(i,:).^senior_par(i,3)+senior_par(i,2)^senior_par(i,3)))*100;
end
SS = SS_dose;
Time = time;
save('senior_neu_ind.mat',"SS",'Time')