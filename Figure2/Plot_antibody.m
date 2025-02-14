clear all
format long
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');

% 95% bootstrap confidence bands
Senior_data = load("senior_ant_ind.mat");
HCW_data = load("hcw_ant_ind.mat");
HCW_fulldose_data = load("hcw_ant_ind_fulldose.mat");
HCW_fulldose = HCW_fulldose_data.H;
HCW = HCW_data.H;
Senior= Senior_data.H;
time = HCW_data.time;
HA_fulldose = squeeze(HCW_fulldose(9,:,:));
HA = squeeze(HCW(9,:,:));
SA = squeeze(Senior(9,:,:));
senior_median = median(SA');
HCW_median = median(HA');
HCW_median_fulldose = median(HA_fulldose');
HH_ci = zeros(2, 4000); % Upper and lower CI bounds for HA data
HH_fulldose_ci = zeros(2, 4000); % Upper and lower CI bounds for HA_fulldose data
SS_ci = zeros(2, 4000); % Upper and lower CI bounds for SA data
for i = 1:4000
    % Extract data for the current simulation
    ha_data = HA(i,:);
    ha_fulldose_data = HA_fulldose(i,:);
    sa_data = SA(i,:);
    
    % Calculate 95% CI using percentiles
    HH_ci(1,i) = prctile(ha_data, 97.5); % Upper bound (97.5th percentile)
    HH_ci(2,i) = prctile(ha_data, 2.5);  % Lower bound (2.5th percentile)
    
    HH_fulldose_ci(1,i) = prctile(ha_fulldose_data, 97.5); % Upper bound (97.5th percentile)
    HH_fulldose_ci(2,i) = prctile(ha_fulldose_data, 2.5);  % Lower bound (2.5th percentile)
    
    SS_ci(1,i) = prctile(sa_data, 97.5); % Upper bound (97.5th percentile)
    SS_ci(2,i) = prctile(sa_data, 2.5);  % Lower bound (2.5th percentile)
end
%%
a1 =max(HCW_median(1001:2000));
a2 =max(senior_median(1001:2000));
a3 =max(HCW_median);
a4 =max(senior_median);
a5 =max(HCW_median_fulldose);
%%
senior_sample = readtable('senior.csv' );
hcw_sample = readtable('hcw.csv' );
% senior_mean = nanmean(table2array(senior_sample(:,2:8)));
% hcw_mean =  nanmean(table2array(hcw_sample(:,2:8)));
% senior_std = nanstd(table2array(senior_sample(:,2:8)));
% hcw_std = nanstd(table2array(hcw_sample(:,2:8)));
senior_median_data = nanmedian(table2array(senior_sample(:,2:8)));
hcw_median_data =  nanmedian(table2array(hcw_sample(:,2:8)));
timepoints = [30,90,150,240,300,360,450];
x = [0 700 700 0];
y1 = [9e4 9e4 2.5e4 2.5e4];
y2 = [0 0 2.5e4 2.5e4];
%
figure
plot(linspace(0, 700, 1e3), repelem(2.5e4, 1000), '--', 'LineWidth', 1.5,'Color', [180,180,180] / 255);
hold on;
patch([time, fliplr(time)], [HH_ci(1, :), fliplr(HH_ci(2, :))], 1, ...
    'FaceColor', '#b24a4b', 'EdgeColor', 'none', 'FaceAlpha', 0.2); % HH CI
hold on;
patch([time, fliplr(time)], [SS_ci(1, :), fliplr(SS_ci(2, :))], 1, ...
    'FaceColor', '#4276bc', 'EdgeColor', 'none', 'FaceAlpha', 0.2); % SS CI
hold on;
scatter(timepoints, senior_median_data, 15, 'filled','MarkerFaceColor','#4276bc') ;
hold on
scatter(timepoints, hcw_median_data, 15, 'filled','MarkerFaceColor','#b24a4b')
hold on
plot(time,senior_median, '-','LineWidth',1.2,'color','#4276bc');
hold on
plot(time,HCW_median,'-','LineWidth',1.2,'color','#b24a4b');
hold on
hold off
box off
xlim([0 630])
ylim([0 7e4])
set(gca,'FontSize',7)
xlabel('Time (Days) ')
ylabel('Antibodies (AU/mL)')
figfile = fullfile(pathname,'Antibody_region');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 7 4.5]); 
set(gcf, 'PaperSize', [7 4.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);


%%

timepoints = [30,90,150,240,300,360,450];
x = [0 700 700 0];
y1 = [1.4e5 1.4e5 2.5e4 2.5e4];
y2 = [0 0 2.5e4 2.5e4];
figure
plot(linspace(0, 700, 1e3), repelem(2.5e4, 1000), '--', 'LineWidth', 1,'Color', [180,180,180] / 255);
hold on;
%plot(linspace(0,700,1e3), repelem(2.5e4,1000),'--','LineWidth',1,'Color',[150,150,150]/255)
hold on
patch([time, fliplr(time)],[HH_fulldose_ci(1,:), fliplr(HH_fulldose_ci(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.2); %CI
hold on
patch([time, fliplr(time)],[SS_ci(1,:), fliplr(SS_ci(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.2); %CI
ylabel('Antibodies (AU/mL)')
hold on
plot(time,senior_median, '-','LineWidth',1.2,'color','#4276bc');
hold on
plot(time,HCW_median_fulldose,'-','LineWidth',1.2,'color','#b24a4b');
box off
xlim([0 630])
ylim([0 1e5])
set(gca,'FontSize',7)
xlabel('Time (Days) ')
figfile = fullfile(pathname,'Antibody_fulldoser_region');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);
%%
data = load('hcw_ant_ind.mat');
H = data.H;
time= data.time;

x = [30, 90, 150,240,300,360,450];

%
hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
anti_hcw = readtable('hcw.csv');
AH = table2array(anti_hcw(:,2:8));
figure(3)
plot(linspace(0,700,1e3), repelem(2.5e4,1000),'--','LineWidth',1,'Color',[200,200,200]/255)
hold on

for i = 1:length(hcw_lam3)
    plot(time,H(9,:,i),'-','LineWidth',0.5,'color','#DDB9B9')
    hold on
end
h22 = scatter (x, AH,4,'o','filled','MarkerFaceColor','#b24a4b') ;
hold off
box off
ylabel('Antibodies (AU/mL)')
ylim([0 7e4])
xlim([0 630])
set(gca,'FontSize',7)
xlabel('Time (Days) ')
figfile = fullfile(pathname,'hcw_antibody_ind');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);


%%
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
data = load('senior_ant_ind.mat');
H = data.H;
time= data.time;
idx30= min(find(time>=30));
idx90= min(find(time>=90));
idx150= min(find(time>=150));
idx240= min(find(time>=240));
idx300= min(find(time>=300));
idx360= min(find(time>=360));
idx450= min(find(time>=450));

x = [30, 90, 150,240,300,360,450];

anti_senior = readtable('senior.csv');
AH = table2array(anti_senior(:,2:8));
figure(4)
plot(linspace(0,700,1e3), repelem(2.5e4,1000),'--','LineWidth',1,'Color',[200,200,200]/255)
hold on
for i = 1:length(senior_dt)
    plot(time,H(9,:,i),'-','LineWidth',0.5,'color','#BDC9DF')
    hold on
end
h2 = scatter (x, AH,4,'o','filled','MarkerFaceColor','#4276bc') ;
hold on
ylabel('Antibodies (AU/mL)')
box off
ylim([0 7e4])
xlim([0 630])
set(gca,'FontSize',7)
xlabel('Time (Days) ')
figfile = fullfile(pathname,'senior_antibody_ind');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%

% 95% bootstrap confidence bands
Senior_data = load("Senior_ant_ind.mat");
HCW_data = load("HCW_ant_ind.mat");
HCW_fulldose_data = load("HCW_ant_ind_fulldose.mat");

time = HCW_data.time;
HCW_fulldose = HCW_fulldose_data.H;
HCW = HCW_data.H;
Senior= Senior_data.H;
HT_fulldose_region = zeros(2,4000);
HT_region = zeros(2,4000);
ST_region = zeros(2,4000);
HB_fulldose_region = zeros(2,4000);
HB_region = zeros(2,4000);
SB_region = zeros(2,4000);
HL_fulldose_region = zeros(2,4000);
HL_region = zeros(2,4000);
SL_region = zeros(2,4000);
HS_fulldose_region = zeros(2,4000);
HS_region = zeros(2,4000);
SS_region = zeros(2,4000);
HM_fulldose_region = zeros(2,4000);
HM_region = zeros(2,4000);
SM_region = zeros(2,4000);




for i = 1:4000
    % T cells
    HT_region(:,i) = prctile(squeeze(HCW(3,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    HT_fulldose_region(:,i) = prctile(squeeze(HCW_fulldose(3,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    ST_region(:,i) = prctile(squeeze(Senior(3,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    %Germinal center B cells
    HB_region(:,i) = prctile(squeeze(HCW(5,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    HB_fulldose_region(:,i) = prctile(squeeze(HCW_fulldose(5,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    SB_region(:,i) = prctile(squeeze(Senior(5,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    %long-lived plasma cells
    HL_region(:,i) = prctile(squeeze(HCW(6,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    HL_fulldose_region(:,i) = prctile(squeeze(HCW_fulldose(6,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    SL_region(:,i) = prctile(squeeze(Senior(6,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
   %short-lived plasma cells
    HS_region(:,i) = prctile(squeeze(HCW(7,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    HS_fulldose_region(:,i) = prctile(squeeze(HCW_fulldose(7,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    SS_region(:,i) = prctile(squeeze(Senior(7,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
   %memory cells
    HM_region(:,i) = prctile(squeeze(HCW(8,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    HM_fulldose_region(:,i) = prctile(squeeze(HCW_fulldose(8,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
    SM_region(:,i) = prctile(squeeze(Senior(8,i,:)),[2.5, 97.5])*1e9; % Calculates the 95% CB at every simulated point
 

end
%%
senior_median = median(Senior,3);
HCW_median = median(HCW,3);
HCW_median_fulldose = median(HCW_fulldose,3);


%%
% full dose
figure
%plot(linspace(0,700,1e3), repelem(2.5e4,1000),'--','LineWidth',1,'Color',[150,150,150]/255)
hold on
patch([time, fliplr(time)],[HT_fulldose_region(1,:), fliplr(HT_fulldose_region(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.2); %CI
hold on
patch([time, fliplr(time)],[ST_region(1,:), fliplr(ST_region(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.2); %CI
hold on
plot(time,senior_median(3,:)*1e9,'-','LineWidth',1.2,'color','#4276bc');
hold on
plot(time,HCW_median_fulldose(3,:)*1e9,'-','LineWidth',1.2,'color','#b24a4b');
hold on
box off
xlim([0 630])
ylim([0 1.6e8])
yticks([0,0.4,0.8,1.2,1.6]*1e8)
set(gca,'FontSize',7)
xlabel('Time (Days) ')
ylabel({'Tfh cells (cells/mL)'})
figfile = fullfile(pathname,'Th_fulldose_region');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);
%%
figure
patch([time, fliplr(time)],[HB_fulldose_region(1,:), fliplr(HB_fulldose_region(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.2); %CI
hold on
patch([time, fliplr(time)],[SB_region(1,:), fliplr(SB_region(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.2); %CI
hold on
plot(time,senior_median(5,:)*1e9,'-','LineWidth',1.2,'color','#4276bc');
hold on
plot(time,HCW_median_fulldose(5,:)*1e9,'-','LineWidth',1.2,'color','#b24a4b');
hold on
box off
xlim([0 630])
ylim([0 4e8])
set(gca,'FontSize',7)
xlabel('Time (Days) ')
ylabel({'GC B cells (cells/mL)'})
figfile = fullfile(pathname,'GB_fulldose_region');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);
%% short-lived
figure
%plot(linspace(0,700,1e3), repelem(2.5e4,1000),'--','LineWidth',1,'Color',[150,150,150]/255)
hold on
patch([time, fliplr(time)],[HS_fulldose_region(1,:), fliplr(HS_fulldose_region(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.2); %CI
hold on
patch([time, fliplr(time)],[SS_region(1,:), fliplr(SS_region(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.2); %CI
hold on
plot(time,senior_median(7,:)*1e9,'-','LineWidth',1.2,'color','#4276bc');
hold on
plot(time,HCW_median_fulldose(7,:)*1e9,'-','LineWidth',1.2,'color','#b24a4b');
hold on
box off
xlim([0 630])
ylim([0 1.5e9])
set(gca,'FontSize',7)
xlabel('Time (Days) ')
ylabel({'Plasmablasts (cells/mL)'})
figfile = fullfile(pathname,'SLPC_fulldose_region');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

%% long-lived
figure
%plot(linspace(0,700,1e3), repelem(2.5e4,1000),'--','LineWidth',1,'Color',[150,150,150]/255)
hold on
patch([time, fliplr(time)],[HL_fulldose_region(1,:), fliplr(HL_fulldose_region(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.2); %CI
hold on
patch([time, fliplr(time)],[SL_region(1,:), fliplr(SL_region(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.2); %CI
hold on
plot(time,senior_median(6,:)*1e9,'-','LineWidth',1.2,'color','#4276bc');
hold on
plot(time,HCW_median_fulldose(6,:)*1e9,'-','LineWidth',1.2,'color','#b24a4b');
hold on
box off
xlim([0 630])
ylim([0  4e9])
set(gca,'FontSize',7)
xlabel('Time (Days) ')
ylabel({'Plasma cells (cells/mL)'})
figfile = fullfile(pathname,'LLPC_fulldose_region');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);
%% memory
figure
%plot(linspace(0,700,1e3), repelem(2.5e4,1000),'--','LineWidth',1,'Color',[150,150,150]/255)
hold on
patch([time, fliplr(time)],[HM_fulldose_region(1,:), fliplr(HM_fulldose_region(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.2); %CI
hold on
patch([time, fliplr(time)],[SM_region(1,:), fliplr(SM_region(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.2); %CI
hold on
plot(time,senior_median(8,:)*1e9,'-','LineWidth',1.2,'color','#4276bc');
hold on
plot(time,HCW_median_fulldose(8,:)*1e9,'-','LineWidth',1.2,'color','#b24a4b');
hold on
box off
xlim([0 630])
ylim([0 1.2e10])
yticks([0, 0.4, 0.8, 1.2]*1e10)
set(gca,'FontSize',7)
xlabel('Time (Days) ')
ylabel('Memory B cells (cells/mL)')
figfile = fullfile(pathname,'Memory_fulldose_region');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);
close all

