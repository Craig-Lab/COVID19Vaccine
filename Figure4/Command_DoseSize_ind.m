%% LLL hcw
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

LLL = zeros(10,4000, length(hcw_lam2));


    dose1 = 0.5;
    dose2 = 0.5;
    dose3 = 0.5;
   
 for j = 1:length(hcw_lam2)
            p.lam2 = hcw_lam2(j);
            p.lam3 = hcw_lam3(j);
            p.d_t  = hcw_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            LLL(:,:,j) = sol;
 end
save('HCW_LLL.mat','time','LLL')

%% LHL hcw
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

LHL = zeros(10,4000, length(hcw_lam2));


    dose1 = 0.5;
    dose2 = 1;
    dose3 = 0.5;
   
 for j = 1:length(hcw_lam2)
            p.lam2 = hcw_lam2(j);
            p.lam3 = hcw_lam3(j);
            p.d_t  = hcw_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            LHL(:,:,j) = sol;
 end
save('HCW_LHL.mat','time','LHL')

%% LLH hcw
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

LLH = zeros(10,4000, length(hcw_lam2));


    dose1 = 0.5;
    dose2 = 0.5;
    dose3 = 1;
   
 for j = 1:length(hcw_lam2)
            p.lam2 = hcw_lam2(j);
            p.lam3 = hcw_lam3(j);
            p.d_t  = hcw_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
           LLH(:,:,j) = sol;
 end
save('HCW_LLH.mat','time','LLH')

%% LHH hcw
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
load_parameters_new()

LHH = zeros(10,4000, length(hcw_lam2));

    dose1 = 0.5;
    dose2 = 1;
    dose3 = 1;
   
 for j = 1:length(hcw_lam2)
            p.lam2 = hcw_lam2(j);
            p.lam3 = hcw_lam3(j);
            p.d_t  = hcw_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            LHH(:,:,j) = sol;
 end
save('HCW_LHH.mat','time','LHH')
%% HLL hcw
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

HLL = zeros(10,4000, length(hcw_lam2));


    dose1 = 1;
    dose2 = 0.5;
    dose3 = 0.5;
   
 for j = 1:length(hcw_lam2)
            p.lam2 = hcw_lam2(j);
            p.lam3 = hcw_lam3(j);
            p.d_t  = hcw_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            HLL(:,:,j) = sol;
 end
save('HCW_HLL.mat','time','HLL')

%% HHL hcw
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

HHL = zeros(10,4000, length(hcw_lam2));


    dose1 = 1;
    dose2 = 1;
    dose3 = 0.5;
   
 for j = 1:length(hcw_lam2)
            p.lam2 = hcw_lam2(j);
            p.lam3 = hcw_lam3(j);
            p.d_t  = hcw_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            HHL(:,:,j) = sol;
 end
save('HCW_HHL.mat','time','HHL')

%% HLH hcw
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

HLH = zeros(10,4000, length(hcw_lam2));


    dose1 = 1;
    dose2 = 0.5;
    dose3 = 1;
   
 for j = 1:length(hcw_lam2)
            p.lam2 = hcw_lam2(j);
            p.lam3 = hcw_lam3(j);
            p.d_t  = hcw_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
           HLH(:,:,j) = sol;
 end
save('HCW_HLH.mat','time','HLH')

%% HHH hcw
clear all
format long

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

HHH = zeros(10,4000, length(hcw_lam2));


    dose1 = 1;
    dose2 = 1;
    dose3 = 1;
   
 for j = 1:length(hcw_lam2)
            p.lam2 = hcw_lam2(j);
            p.lam3 = hcw_lam3(j);
            p.d_t  = hcw_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            HHH(:,:,j) = sol;
 end
save('HCW_HHH.mat','time','HHH')

%% LLL senior
clear all
format long

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

LLL = zeros(10,4000, length(senior_lam2));

    dose1 = 0.5;
    dose2 = 0.5;
    dose3 = 0.5;
   
 for j = 1:length(senior_lam2)
            p.lam2 = senior_lam2(j);
            p.lam3 = senior_lam3(j);
            p.d_t  = senior_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            LLL(:,:,j) = sol;
 end
save('senior_LLL.mat','time','LLL')

%% LHL senior
clear all
format long

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

LHL = zeros(10,4000, length(senior_lam2));


    dose1 = 0.5;
    dose2 = 1;
    dose3 = 0.5;
   
 for j = 1:length(senior_lam2)
            p.lam2 = senior_lam2(j);
            p.lam3 = senior_lam3(j);
            p.d_t  = senior_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            LHL(:,:,j) = sol;
 end
save('senior_LHL.mat','time','LHL')

%% LLH senior
clear all
format long

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

LLH = zeros(10,4000, length(senior_lam2));


    dose1 = 0.5;
    dose2 = 0.5;
    dose3 = 1;
   
 for j = 1:length(senior_lam2)
            p.lam2 = senior_lam2(j);
            p.lam3 = senior_lam3(j);
            p.d_t  = senior_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
           LLH(:,:,j) = sol;
 end
save('senior_LLH.mat','time','LLH')

%% LHH senior
clear all
format long

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

LHH = zeros(10,4000, length(senior_lam2));


    dose1 = 0.5;
    dose2 = 1;
    dose3 = 1;
   
 for j = 1:length(senior_lam2)
            p.lam2 = senior_lam2(j);
            p.lam3 = senior_lam3(j);
            p.d_t  = senior_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            LHH(:,:,j) = sol;
 end
save('senior_LHH.mat','time','LHH')

%% HLL senior
clear all
format long

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

HLL = zeros(10,4000, length(senior_lam2));


    dose1 = 1;
    dose2 = 0.5;
    dose3 = 0.5;
   
 for j = 1:length(senior_lam2)
            p.lam2 = senior_lam2(j);
            p.lam3 = senior_lam3(j);
            p.d_t  = senior_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            HLL(:,:,j) = sol;
 end
save('senior_HLL.mat','time','HLL')

%% HHL senior
clear all
format long

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

HHL = zeros(10,4000, length(senior_lam2));


    dose1 = 1;
    dose2 = 1;
    dose3 = 0.5;
   
 for j = 1:length(senior_lam2)
            p.lam2 = senior_lam2(j);
            p.lam3 = senior_lam3(j);
            p.d_t  = senior_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            HHL(:,:,j) = sol;
 end
save('senior_HHL.mat','time','HHL')

%% HLH senior
clear all
format long

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

HLH = zeros(10,4000, length(senior_lam2));


    dose1 = 1;
    dose2 = 0.5;
    dose3 = 1;
   
 for j = 1:length(senior_lam2)
            p.lam2 = senior_lam2(j);
            p.lam3 = senior_lam3(j);
            p.d_t  = senior_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
           HLH(:,:,j) = sol;
 end
save('senior_HLH.mat','time','HLH')

%% HHH senior
clear all
format long

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
%
% lipid nanoparticles
load_parameters_new()

HHH = zeros(10,4000, length(senior_lam2));


    dose1 = 1;
    dose2 = 1;
    dose3 = 1;
   
 for j = 1:length(senior_lam2)
            p.lam2 = senior_lam2(j);
            p.lam3 = senior_lam3(j);
            p.d_t  = senior_dt(j);

            tic
             [sol, time] = Model_3doses_amount(p,dose1, dose2, dose3);
            toc
            HHH(:,:,j) = sol;
 end
save('senior_HHH.mat','time','HHH')
%%
%
hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
load('HCW_LLL.mat')
LLL_neu = zeros(length(hcw_lam2),4000);
hcw_par_data = readtable('hcw_neu_par.xlsx');
par_h = table2array(hcw_par_data(:,2:4));
DS = squeeze(LLL(9,:,:))'./1e3;

for j = 1:length(hcw_lam2)
   LLL_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_HCW_LLL.mat','LLL_neu','time')
%
load('HCW_LHL.mat')
LHL_neu = zeros(length(hcw_lam2),4000);
hcw_par_data = readtable('hcw_neu_par.xlsx');
par_h = table2array(hcw_par_data(:,2:4));
DS = squeeze(LHL(9,:,:))'./1e3;

for j = 1:length(hcw_lam2)
   LHL_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_HCW_LHL.mat','LHL_neu','time')
%
load('HCW_LLH.mat')
LLH_neu = zeros(length(hcw_lam2),4000);
hcw_par_data = readtable('hcw_neu_par.xlsx');
par_h = table2array(hcw_par_data(:,2:4));
DS = squeeze(LLH(9,:,:))'./1e3;
for j = 1:length(hcw_lam2)
   LLH_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_HCW_LLH.mat','LLH_neu','time')
%
load('HCW_LHH.mat')
LHH_neu = zeros(length(hcw_lam2),4000);
hcw_par_data = readtable('hcw_neu_par.xlsx');
par_h = table2array(hcw_par_data(:,2:4));
DS = squeeze(LHH(9,:,:))'./1e3;
for j = 1:length(hcw_lam2)
   LHH_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_HCW_LHH.mat','LHH_neu','time')

%
load('HCW_HLL.mat')
HLL_neu = zeros(length(hcw_lam2),4000);
hcw_par_data = readtable('hcw_neu_par.xlsx');
par_h = table2array(hcw_par_data(:,2:4));
DS = squeeze(HLL(9,:,:))'./1e3;
for j = 1:length(hcw_lam2)
   HLL_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_HCW_HLL.mat','HLL_neu','time')

load('HCW_HHL.mat')
HHL_neu = zeros(length(hcw_lam2),4000);
hcw_par_data = readtable('hcw_neu_par.xlsx');
par_h = table2array(hcw_par_data(:,2:4));
DS = squeeze(HHL(9,:,:))'./1e3;

for j = 1:length(hcw_lam2)
   HHL_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_HCW_HHL.mat','HHL_neu','time')
%
load('HCW_HLH.mat')
HLH_neu = zeros(length(hcw_lam2),4000);
hcw_par_data = readtable('hcw_neu_par.xlsx');
par_h = table2array(hcw_par_data(:,2:4));
DS = squeeze(HLH(9,:,:))'./1e3;

for j = 1:length(hcw_lam2)
   HLH_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_HCW_HLH.mat','HLH_neu','time')

%
load('HCW_HHH.mat')
HHH_neu = zeros(length(hcw_lam2),4000);
hcw_par_data = readtable('hcw_neu_par.xlsx');
par_h = table2array(hcw_par_data(:,2:4));
DS = squeeze(HHH(9,:,:))'./1e3;
for j = 1:length(hcw_lam2)
   HHH_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_HCW_HHH.mat','HHH_neu','time')

%
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
load('senior_LLL.mat')
LLL_neu = zeros(length(senior_lam2),4000);
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
DS = squeeze(LLL(9,:,:))'./1e3;

for j = 1:length(senior_lam2)
   LLL_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_senior_LLL.mat','LLL_neu','time')
%
load('senior_LHL.mat')
LHL_neu = zeros(length(senior_lam2),4000);
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
DS = squeeze(LHL(9,:,:))'./1e3;

for j = 1:length(senior_lam2)
   LHL_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_senior_LHL.mat','LHL_neu','time')
%
load('senior_LLH.mat')
LLH_neu = zeros(length(senior_lam2),4000);
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
DS = squeeze(LLH(9,:,:))'./1e3;
for j = 1:length(senior_lam2)
   LLH_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_senior_LLH.mat','LLH_neu','time')
%
load('senior_LHH.mat')
LHH_neu = zeros(length(senior_lam2),4000);
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
DS = squeeze(LHH(9,:,:))'./1e3;
for j = 1:length(senior_lam2)
   LHH_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_senior_LHH.mat','LHH_neu','time')

%
load('senior_HLL.mat')
HLL_neu = zeros(length(senior_lam2),4000);
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
DS = squeeze(HLL(9,:,:))'./1e3;
for j = 1:length(senior_lam2)
   HLL_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_senior_HLL.mat','HLL_neu','time')

load('senior_HHL.mat')
HHL_neu = zeros(length(senior_lam2),4000);
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
DS = squeeze(HHL(9,:,:))'./1e3;

for j = 1:length(senior_lam2)
   HHL_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_senior_HHL.mat','HHL_neu','time')
%
load('senior_HLH.mat')
HLH_neu = zeros(length(senior_lam2),4000);
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
DS = squeeze(HLH(9,:,:))'./1e3;

for j = 1:length(senior_lam2)
   HLH_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_senior_HLH.mat','HLH_neu','time')

%
load('senior_HHH.mat')
HHH_neu = zeros(length(senior_lam2),4000);
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
DS = squeeze(HHH(9,:,:))'./1e3;
for j = 1:length(senior_lam2)
   HHH_neu(j,:) =real(par_h(j,1)+(1-par_h(j,1))* DS(j,:).^par_h(j,3)./(DS(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('neu_senior_HHH.mat','HHH_neu','time')
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
