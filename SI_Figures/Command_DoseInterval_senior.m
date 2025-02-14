%%
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
clear all
format long
load_parameters_new()
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));

dose3 = 1;
SA = zeros(10,4000,length(senior_lam2));
month1 = 1;
for i = 1:length(senior_lam2)
        p.lam2 = senior_lam2(i);
        p.lam3=senior_lam3(i);
        p.d_t = senior_dt(i);
        delay1 = month1*30;
        [sol, time] = Model_3doses_interval(p,dose3,delay1, 210);
        SA(:,:,i) = sol;
        i
end
save('Dose_1m_senior.mat','time',"SA")

senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
LA = squeeze(SA(9,:,:))'./1e3;
SH = zeros(length(senior_lam2),4000);
for j = 1:length(senior_lam2)
   SH(j,:) =real(par_h(j,1)+(1-par_h(j,1))* LA(j,:).^par_h(j,3)./(LA(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('Neu_senior_int12_1m.mat','SH','time')

%%
clear all
format long
load_parameters_new()
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));

dose3 = 1;
SA = zeros(10,4000,length(senior_lam2));
month1 = 2;
for i = 1:length(senior_lam2)
        p.lam2 = senior_lam2(i);
        p.lam3=senior_lam3(i);
        p.d_t = senior_dt(i);
        delay1 = month1*30;
        [sol, time] = Model_3doses_interval(p,dose3,delay1, 210);
        SA(:,:,i) = sol;
        i
end
save('Dose_2m_senior.mat','time',"SA")

senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
LA = squeeze(SA(9,:,:))'./1e3;
SH = zeros(length(senior_lam2),4000);
for j = 1:length(senior_lam2)
   SH(j,:) =real(par_h(j,1)+(1-par_h(j,1))* LA(j,:).^par_h(j,3)./(LA(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('Neu_senior_int12_2m.mat','SH','time')
clear all
format long
load_parameters_new()
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));

dose3 = 1;
SA = zeros(10,4000,length(senior_lam2));
month1 = 4;
for i = 1:length(senior_lam2)
        p.lam2 = senior_lam2(i);
        p.lam3=senior_lam3(i);
        p.d_t = senior_dt(i);
        delay1 = month1*30;
        [sol, time] = Model_3doses_interval(p,dose3,delay1, 210);
        SA(:,:,i) = sol;
        i
end
save('Dose_4m_senior.mat','time',"SA")

senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
LA = squeeze(SA(9,:,:))'./1e3;
SH = zeros(length(senior_lam2),4000);
for j = 1:length(senior_lam2)
   SH(j,:) =real(par_h(j,1)+(1-par_h(j,1))* LA(j,:).^par_h(j,3)./(LA(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('Neu_senior_int12_4m.mat','SH','time')
%%

clear all
format long
load_parameters_new()
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));

dose3 = 1;
SA = zeros(10,4000,length(senior_lam2));
month1 = 6;
for i = 1:length(senior_lam2)
        p.lam2 = senior_lam2(i);
        p.lam3=senior_lam3(i);
        p.d_t = senior_dt(i);
        delay1 = month1*30;
        [sol, time] = Model_3doses_interval(p,dose3,delay1, 210);
        SA(:,:,i) = sol;
        i
end
save('Dose_6m_senior.mat','time',"SA")

senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
LA = squeeze(SA(9,:,:))'./1e3;
SH = zeros(length(senior_lam2),4000);
for j = 1:length(senior_lam2)
   SH(j,:) =real(par_h(j,1)+(1-par_h(j,1))* LA(j,:).^par_h(j,3)./(LA(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('Neu_senior_int12_6m.mat','SH','time')

pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
clear all
format long
load_parameters_new()
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));

dose3 = 1;
SA = zeros(10,4000,length(senior_lam2));
month2 = 6;
for i = 1:length(senior_lam2)
        p.lam2 = senior_lam2(i);
        p.lam3=senior_lam3(i);
        p.d_t = senior_dt(i);
        delay2 = month2*30;
        [sol, time] = Model_3doses_interval(p,dose3,60,delay2);
        SA(:,:,i) = sol;
        i
end
save('Dose23_6m_senior.mat','time',"SA")

senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
LA = squeeze(SA(9,:,:))'./1e3;
SH = zeros(length(senior_lam2),4000);
for j = 1:length(senior_lam2)
   SH(j,:) =real(par_h(j,1)+(1-par_h(j,1))* LA(j,:).^par_h(j,3)./(LA(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('Neu_senior_int23_6m.mat','SH','time')

%%
clear all
format long
load_parameters_new()
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));

dose3 = 1;
SA = zeros(10,4000,length(senior_lam2));
month2 = 7;
for i = 1:length(senior_lam2)
        p.lam2 = senior_lam2(i);
        p.lam3=senior_lam3(i);
        p.d_t = senior_dt(i);
        delay2 = month2*30;
        [sol, time] = Model_3doses_interval(p,dose3,60,delay2);
        SA(:,:,i) = sol;
        i
end
save('Dose23_7m_senior.mat','time',"SA")

senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
LA = squeeze(SA(9,:,:))'./1e3;
SH = zeros(length(senior_lam2),4000);
for j = 1:length(senior_lam2)
   SH(j,:) =real(par_h(j,1)+(1-par_h(j,1))* LA(j,:).^par_h(j,3)./(LA(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('Neu_senior_int23_7m.mat','SH','time')
%%
clear all
format long
load_parameters_new()
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));

dose3 = 1;
SA = zeros(10,4000,length(senior_lam2));
month2 = 8;
for i = 1:length(senior_lam2)
        p.lam2 = senior_lam2(i);
        p.lam3=senior_lam3(i);
        p.d_t = senior_dt(i);
        delay2 = month2*30;
        [sol, time] = Model_3doses_interval(p,dose3,60,delay2);
        SA(:,:,i) = sol;
        i
end
save('Dose23_8m_senior.mat','time',"SA")

senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
LA = squeeze(SA(9,:,:))'./1e3;
SH = zeros(length(senior_lam2),4000);
for j = 1:length(senior_lam2)
   SH(j,:) =real(par_h(j,1)+(1-par_h(j,1))* LA(j,:).^par_h(j,3)./(LA(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('Neu_senior_int23_8m.mat','SH','time')
%%

clear all
format long
load_parameters_new()
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));

dose3 = 1;
SA = zeros(10,4000,length(senior_lam2));
month2 = 9;
for i = 1:length(senior_lam2)
        p.lam2 = senior_lam2(i);
        p.lam3=senior_lam3(i);
        p.d_t = senior_dt(i);
        delay2 = month2*30;
        [sol, time] = Model_3doses_interval(p,dose3,60,delay2);
        SA(:,:,i) = sol;
        i
end
save('Dose23_9m_senior.mat','time',"SA")

senior_par_data = readtable('senior_neu_par.xlsx');
par_h = table2array(senior_par_data(:,2:4));
LA = squeeze(SA(9,:,:))'./1e3;
SH = zeros(length(senior_lam2),4000);
for j = 1:length(senior_lam2)
   SH(j,:) =real(par_h(j,1)+(1-par_h(j,1))* LA(j,:).^par_h(j,3)./(LA(j,:).^par_h(j,3)+par_h(j,2)^par_h(j,3)))*100;
end
save('Neu_senior_int23_9m.mat','SH','time')






