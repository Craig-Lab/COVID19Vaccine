%%
clear all
load('senior_booster05_yearly.mat')
load('senior_booster1_yearly.mat')
load('hcw_booster05_yearly.mat')
load('hcw_booster1_yearly.mat')
load('senior_booster05_6m.mat')
hcw_half = load('HCW_average_repeated_halfdose.mat').hcw_average_sol;
hcw_full = load('HCW_average_repeated_fulldose.mat').hcw_average_sol;
time_booster = load('HCW_average_repeated_fulldose.mat').hcw_average_time;
senior_half = load('senior_average_repeated_halfdose.mat').senior_average_sol;
senior_full = load('senior_average_repeated_fulldose.mat').senior_average_sol;
senior_half_6m = load('senior_average_repeated_halfdose_6m.mat').senior_average_sol;
load('average_neu_par_monolix.mat')
hcw = par_monolix(1,:);
senior = par_monolix(2,:);
senior_ant_full = senior_full(9,:)./1e3;
senior_ant_half = senior_half(9,:)./1e3;
senior_ant_half_6m = senior_half_6m(9,:)./1e3; 
hcw_ant_full = hcw_full(9,:)./1e3;
hcw_ant_half = hcw_half(9,:)./1e3;
neu_senior_full = (senior (1)+(1-senior(1))*(senior_ant_full .^senior(3))./(senior(2).^senior(3)+senior_ant_full .^senior(3)))*100;
neu_senior_half = (senior (1)+(1-senior(1))*(senior_ant_half.^senior(3))./(senior(2).^senior(3)+senior_ant_half.^senior(3)))*100;
neu_senior_half_6m = (senior (1)+(1-senior(1))*(senior_ant_half_6m.^senior(3))./(senior(2).^senior(3)+senior_ant_half_6m.^senior(3)))*100;

neu_hcw_full = (hcw (1)+(1-hcw(1))*(hcw_ant_full.^hcw(3))./(hcw(2).^hcw(3)+hcw_ant_full.^hcw(3)))*100;
neu_hcw_half = (hcw (1)+(1-hcw(1))*(hcw_ant_half.^hcw(3))./(hcw(2).^hcw(3)+hcw_ant_half.^hcw(3)))*100;

hcw_par_data = readtable('hcw_neu_par.xlsx');
hcw_par = table2array(hcw_par_data(:,2:4));
LH1_data =squeeze(hcw_booster1(9,:,:))'./1e3;
LH05_data =squeeze(hcw_booster05(9,:,:))'./1e3;
LH1 = LH1_data;
LH05 = LH05_data;
HH_neu_booster1 = zeros(length(hcw_par),length(time_booster));
HH_neu_booster05 = zeros(length(hcw_par),length(time_booster));

for i = 1:32
    HH_neu_booster1(i,:) =100*real(hcw_par(i,1)+(1-hcw_par(i,1))* LH1(i,:).^hcw_par(i,3)./(LH1(i,:).^hcw_par(i,3)+hcw_par(i,2)^hcw_par(i,3)));
    HH_neu_booster05(i,:) =100*real(hcw_par(i,1)+(1-hcw_par(i,1))* LH05(i,:).^hcw_par(i,3)./(LH05(i,:).^hcw_par(i,3)+hcw_par(i,2)^hcw_par(i,3)));
end

senior_par_data = readtable('senior_neu_par.xlsx');
senior_par = table2array(senior_par_data(:,2:4));
LS1_data =squeeze(senior_booster1(9,:,:))'./1e3;
LS05_data =squeeze(senior_booster05(9,:,:))'./1e3;
LS05_6m_data =squeeze(senior_booster05_6m(9,:,:))'./1e3;
LS1 = LS1_data;
LS05 = LS05_data;
LS05_6m = LS05_6m_data;
SS_neu_booster1 = zeros(length(senior_par),length(time_booster));
SS_neu_booster05 = zeros(length(senior_par),length(time_booster));
SS_neu_booster05_6m = zeros(length(senior_par),length(time_booster05_6m));
for i = 1:27
    SS_neu_booster1(i,:) =100*real(senior_par(i,1)+(1-senior_par(i,1))* LS1(i,:).^senior_par(i,3)./(LS1(i,:).^senior_par(i,3)+senior_par(i,2)^senior_par(i,3)));
    SS_neu_booster05(i,:) =100*real(senior_par(i,1)+(1-senior_par(i,1))* LS05(i,:).^senior_par(i,3)./(LS05(i,:).^senior_par(i,3)+senior_par(i,2)^senior_par(i,3)));
    SS_neu_booster05_6m(i,:) =100*real(senior_par(i,1)+(1-senior_par(i,1))* LS05_6m(i,:).^senior_par(i,3)./(LS05_6m(i,:).^senior_par(i,3)+senior_par(i,2)^senior_par(i,3)));

end

HH_ant_booster1_ci = zeros(2,length(time_booster));
SS_ant_booster1_ci = zeros(2,length(time_booster));
HH_ant_booster1 = squeeze(hcw_booster1(9,:,:))';
SS_ant_booster1 = squeeze(senior_booster1(9,:,:))';

HH_ant_booster05_ci = zeros(2,length(time_booster));
SS_ant_booster05_ci = zeros(2,length(time_booster));
HH_ant_booster05 = squeeze(hcw_booster05(9,:,:))';
SS_ant_booster05 = squeeze(senior_booster05(9,:,:))';

HH_neu_booster1_ci = zeros(2,length(time_booster));
SS_neu_booster1_ci = zeros(2,length(time_booster));

HH_neu_booster05_ci = zeros(2,length(time_booster));
SS_neu_booster05_ci = zeros(2,length(time_booster));
SS_neu_booster05_6m_ci = zeros(2,length(time_booster05_6m));


for i = 1:length(time_booster05_6m)
      SS_neu_booster05_6m_ci(:,i) = prctile(SS_neu_booster05_6m(:,i),[2.5, 97.5]); % Calculates the 95% CB at every simulated point    
end
for i = 1:length(time_booster)
        HH_neu_booster1_ci(:,i) = prctile(HH_neu_booster1(:,i) ,[2.5, 97.5]); % Calculates the 95% CB at every simulated point
        SS_neu_booster1_ci(:,i) = prctile(SS_neu_booster1(:,i),[2.5, 97.5]); % Calculates the 95% CB at every simulated point
        HH_neu_booster05_ci(:,i) = prctile(HH_neu_booster05(:,i) ,[2.5, 97.5]); % Calculates the 95% CB at every simulated point
        SS_neu_booster05_ci(:,i) = prctile(SS_neu_booster05(:,i),[2.5, 97.5]); % Calculates the 95% CB at every simulated point
        HH_ant_booster1_ci(:,i) = prctile(HH_ant_booster1(:,i) ,[2.5, 97.5]); % Calculates the 95% CB at every simulated point
        SS_ant_booster1_ci(:,i) = prctile(SS_ant_booster1(:,i),[2.5, 97.5]); % Calculates the 95% CB at every simulated point
        HH_ant_booster05_ci(:,i) = prctile(HH_ant_booster05(:,i) ,[2.5, 97.5]); % Calculates the 95% CB at every simulated point
        SS_ant_booster05_ci(:,i) = prctile(SS_ant_booster05(:,i),[2.5, 97.5]); % Calculates the 95% CB at every simulated point
end
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
%%
map= {'#F8F9FA','#E9ECEF','#DEE2E6','#CED4DA','#ADB5BD'};
col_days = validatecolor(map, 'multiple');
x = [0 time_booster(8000) time_booster(8000) 0];
y1 =[70 70 60 60] ;
y2 = [80 80 70 70];
y3 =  [90 90 80 80];
y4 = [95 95 90 90];
y5 = [100 100 95 95];
figure(2)
patch(x,y1,col_days(1,:),'LineStyle','none')
hold on
patch(x,y2,col_days(2,:),'LineStyle','none')
hold on
patch(x,y3,col_days(3,:),'LineStyle','none')
hold on
patch(x,y4,col_days(4,:),'LineStyle','none')
hold on
patch(x,y5,col_days(5,:),'LineStyle','none')
hold on
patch([time_booster, fliplr(time_booster)],[SS_neu_booster1_ci(1,:), fliplr(SS_neu_booster1_ci(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.15); %ci
hold on
patch([time_booster, fliplr(time_booster)],[HH_neu_booster1_ci(1,:), fliplr(HH_neu_booster1_ci(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.15); %ci
hold on
%plot(time_booster,neu_senior_full,'-','LineWidth',1.5,'color','#4276bc');
plot(time_booster,median(SS_neu_booster1),'-','LineWidth',1.2,'color','#4276bc');
hold on
%plot(time_booster,neu_hcw_ful,'-','LineWidth',1.5,'color','#4276bc');
plot(time_booster,median(HH_neu_booster1),'-','LineWidth',1.2,'color','#b24a4b');
hold off
%plot(xlim, [60 60], '--', 'LineWidth', 1.5, 'Color', 'k');
%hold on
xlim([0 time_booster(8000)])
ylim([0 100])
xticks([60,270,630,990,1350,1710,2070,2430,2790])
xlabel('Time (Days) ')
ylabel('Neutralization (%)')
box off
figfile = fullfile(pathname,'Neutralization_booster_full');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 7.5 3.5]); 
set(gcf, 'PaperSize', [7.5 3.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%
figure(4)
patch(x,y1,col_days(1,:),'LineStyle','none')
hold on
patch(x,y2,col_days(2,:),'LineStyle','none')
hold on
patch(x,y3,col_days(3,:),'LineStyle','none')
hold on
patch(x,y4,col_days(4,:),'LineStyle','none')
hold on
patch(x,y5,col_days(5,:),'LineStyle','none')
hold on
patch([time_booster, fliplr(time_booster)],[SS_neu_booster05_ci(1,:), fliplr(SS_neu_booster05_ci(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.15); %ci
hold on
patch([time_booster, fliplr(time_booster)],[HH_neu_booster05_ci(1,:), fliplr(HH_neu_booster05_ci(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.15); %ci
hold on
%plot(time_booster,neu_senior_half,'-','LineWidth',1.5,'color','#4276bc');
plot(time_booster,median(SS_neu_booster05),'-','LineWidth',1.2,'color','#4276bc');
hold on
%plot(time_booster,neu_senior_half,'-','LineWidth',1.5,'color','#4276bc');
plot(time_booster,median(HH_neu_booster05),'-','LineWidth',1.2,'color','#b24a4b');
hold on
xlim([0 time_booster(8000)])
ylim([0 100])
xticks([60,270,630,990,1350,1710,2070,2430,2790])
xlabel('Time (Days) ')
ylabel('Neutralization (%)')
box off
figfile = fullfile(pathname,'Neutralization_booster_half');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 7.5 3.5]); 
set(gcf, 'PaperSize', [7.5 3.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);


%
figure(6)
patch(x,y1,col_days(1,:),'LineStyle','none')
hold on
patch(x,y2,col_days(2,:),'LineStyle','none')
hold on
patch(x,y3,col_days(3,:),'LineStyle','none')
hold on
patch(x,y4,col_days(4,:),'LineStyle','none')
hold on
patch(x,y5,col_days(5,:),'LineStyle','none')
hold on
patch([time_booster, fliplr(time_booster)],[SS_neu_booster1_ci(1,:), fliplr(SS_neu_booster1_ci(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.15); %ci
hold on
patch([time_booster, fliplr(time_booster)],[HH_neu_booster05_ci(1,:), fliplr(HH_neu_booster05_ci(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.15); %ci
hold on
%plot(time_booster,neu_senior_full,'-','LineWidth',1.5,'color','#4276bc');
plot(time_booster,median(SS_neu_booster1),'-','LineWidth',1.2,'color','#4276bc');
hold on
%plot(time_booster, neu_hcw_half,'-','LineWidth',1.5,'color','#b24a4b');
plot(time_booster, median(HH_neu_booster05),'-','LineWidth',1.2,'color','#b24a4b');
hold on
xlim([0 time_booster(8000)])
ylim([0 100])
xticks([60,270,630,990,1350,1710,2070,2430,2790])
xlabel('Time (Days) ')
ylabel('Neutralization (%)')
box off
figfile = fullfile(pathname,'Neutralization_booster_fullhalf');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 7.5 4]); 
set(gcf, 'PaperSize', [7.5 4]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%
map= {'#F8F9FA','#E9ECEF','#DEE2E6','#CED4DA','#ADB5BD'};
col_days = validatecolor(map, 'multiple');
x = [0 time_booster(8000) time_booster(8000) 0];
%%
figure(8)
patch(x,y1,col_days(1,:),'LineStyle','none')
hold on
patch(x,y2,col_days(2,:),'LineStyle','none')
hold on
patch(x,y3,col_days(3,:),'LineStyle','none')
hold on
patch(x,y4,col_days(4,:),'LineStyle','none')
hold on
patch(x,y5,col_days(5,:),'LineStyle','none')
hold on
patch([time_booster, fliplr(time_booster)],[SS_neu_booster1_ci(1,:), fliplr(SS_neu_booster1_ci(2,:))],1,'facecolor','#4276bc','EdgeAlpha',0,'facealpha', 0.15);
hold on
patch([time_booster05_6m, fliplr(time_booster05_6m)],[SS_neu_booster05_6m_ci(1,:), fliplr(SS_neu_booster05_6m_ci(2,:))],1,'FaceColor','#3E7748','EdgeAlpha',0,'facealpha', 0.15);
hold on
%plot(time_booster,neu_senior_full,'-','LineWidth',1.5,'color','#4276bc');
plot(time_booster,median(SS_neu_booster1),'-','LineWidth',1.2,'color','#4276bc');
hold on
plot(time_booster05_6m, median(SS_neu_booster05_6m),'-','LineWidth',1.2,'color','#3E7748');
hold on
xlim([0 time_booster(8000)])
ylim([0 100])
xticks([60,270,630,990,1350,1710,2070,2430,2790])
xlabel('Time (Days) ')
ylabel('Neutralization (%)')
box off

% Set the figure and paper position to increase resolution
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 7.5 4]); 
set(gcf, 'PaperSize', [7.5 4]);
set(gca, 'LooseInset', get(gca,'TightInset'))

% Save the figure
figfile = fullfile(pathname,'Neutralization_booster_1y6m');
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);



idx = find(time_booster05_6m>1350+360,1);
a = median(SS_neu_booster05_6m);
a(idx)
%%
Edges = 0:10:100;
figure(51)
h1 = histogram(HH_neu_booster05(:,5000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster05(:,5000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
%histfit(HH_neu_booster05(:,5000))
hold on
%histfit(SS_neu_booster05(:,5000))
hold on
plot([median(HH_neu_booster05(:,5000)), median(HH_neu_booster05(:,5000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster05(:,5000)), median(SS_neu_booster05(:,5000))],[0 h2.Values(7)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_4_half');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%%
figure(52)
h1 = histogram(HH_neu_booster05(:,6000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster05(:,6000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,6000)), median(HH_neu_booster05(:,6000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster05(:,6000)), median(SS_neu_booster05(:,6000))],[0 h2.Values(6)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_5_half');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%
figure(53)
h1 = histogram(HH_neu_booster05(:,7000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster05(:,7000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,7000)), median(HH_neu_booster05(:,7000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster05(:,7000)), median(SS_neu_booster05(:,7000))],[0 h2.Values(6)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_5_half');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%
figure(54)
h1 = histogram(HH_neu_booster05(:,8000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster05(:,8000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,8000)), median(HH_neu_booster05(:,8000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster05(:,8000)), median(SS_neu_booster05(:,8000))],[0 h2.Values(6)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_6_half');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%
figure(55)
h1 = histogram(HH_neu_booster05(:,9000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster05(:,9000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,9000)), median(HH_neu_booster05(:,9000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster05(:,9000)), median(SS_neu_booster05(:,9000))],[0 h2.Values(6)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_7_half');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);


figure(56)
h1 = histogram(HH_neu_booster05(:,10000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster05(:,10000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,10000)), median(HH_neu_booster05(:,10000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster05(:,10000)), median(SS_neu_booster05(:,10000))],[0 h2.Values(6)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_8_half');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);


close all
%%
Edges = 0:10:100;
figure(61)
h1 = histogram(HH_neu_booster05(:,5000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,5000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,5000)), median(HH_neu_booster05(:,5000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,5000)), median(SS_neu_booster1(:,5000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_4_fullhalf');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%
figure(62)
h1 = histogram(HH_neu_booster05(:,6000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,6000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,6000)), median(HH_neu_booster05(:,6000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,6000)), median(SS_neu_booster1(:,6000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_5_fullhalf');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%
figure(63)
h1 = histogram(HH_neu_booster05(:,7000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,7000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,7000)), median(HH_neu_booster05(:,7000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,7000)), median(SS_neu_booster1(:,7000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_5_fullhalf');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

figure(64)
h1 = histogram(HH_neu_booster05(:,8000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,8000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,8000)), median(HH_neu_booster05(:,8000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,8000)), median(SS_neu_booster1(:,8000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_6_fullhalf');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%
figure(65)
h1 = histogram(HH_neu_booster05(:,9000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,9000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,9000)), median(HH_neu_booster05(:,9000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,9000)), median(SS_neu_booster1(:,9000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_7_fullhalf');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

Edges = 0:10:100;
figure(66)
h1 = histogram(HH_neu_booster05(:,10000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,10000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster05(:,10000)), median(HH_neu_booster05(:,10000))],[0 h1.Values(8)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,10000)), median(SS_neu_booster1(:,10000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_8_fullhalf');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%%
Edges = 0:10:100;
figure(71)
h1 = histogram(HH_neu_booster1(:,5000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,5000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster1(:,5000)), median(HH_neu_booster1(:,5000))],[0 h1.Values(9)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,5000)), median(SS_neu_booster1(:,5000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.5])
xlim([20 100])
yticks([0 0.2 0.5 ])
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_4_full');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%%
figure(72)
h1 = histogram(HH_neu_booster1(:,6000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,6000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster1(:,6000)), median(HH_neu_booster1(:,6000))],[0 h1.Values(9)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,6000)), median(SS_neu_booster1(:,6000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.5])
xlim([20 100])
yticks([0 0.2 0.5])
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_5_full');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%
figure(73)
h1 = histogram(HH_neu_booster1(:,7000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,7000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster1(:,7000)), median(HH_neu_booster1(:,7000))],[0 h1.Values(9)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,7000)), median(SS_neu_booster1(:,7000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.5])
xlim([20 100])
yticks([0 0.2 0.5])
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_5_full');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%%
figure(74)
h1 = histogram(HH_neu_booster1(:,8000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,8000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster1(:,8000)), median(HH_neu_booster1(:,8000))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,8000)), median(SS_neu_booster1(:,8000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.5])
xlim([20 100])
yticks([0 0.2 0.5])
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_6_full');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%
figure(75)
h1 = histogram(HH_neu_booster1(:,9000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,9000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster1(:,9000)), median(HH_neu_booster1(:,9000))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,9000)), median(SS_neu_booster1(:,9000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.5])
xlim([20 100])
yticks([0 0.2 0.5])
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_7_full');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%
Edges = 0:10:100;
figure(76)
h1 = histogram(HH_neu_booster1(:,10000),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,10000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(HH_neu_booster1(:,10000)), median(HH_neu_booster1(:,10000))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(SS_neu_booster1(:,10000)), median(SS_neu_booster1(:,10000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.5])
xlim([20 100])
yticks([0 0.2  0.5])
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_8_full');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%
Edges = 0:10:100;
figure(81)
h1 = histogram(SS_neu_booster05_6m(:,6000),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.3);
hold on %#3e7749');%'#3e7749'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,5000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(SS_neu_booster05_6m(:,6000)), median(SS_neu_booster05_6m(:,6000))],[0 h1.Values(9)],'LineWidth',1.5,'Color','#3e7749');
hold on
plot([median(SS_neu_booster1(:,5000)), median(SS_neu_booster1(:,5000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
ylim([0 0.6])
xlim([20 100])
yticks(0:0.2:0.6)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_4_1y6m');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%
figure(82)
h1 = histogram(SS_neu_booster05_6m(:,8000),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.3);
hold on %#3e7749');%'#3e7749'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,6000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(SS_neu_booster05_6m(:,8000)), median(SS_neu_booster05_6m(:,8000))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#3e7749');
hold on
plot([median(SS_neu_booster1(:,6000)), median(SS_neu_booster1(:,6000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');

xlim([20 100])
yticks(0:0.2:0.6)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_5_1y6m');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%
figure(83)
h1 = histogram(SS_neu_booster05_6m(:,10000),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.3);
hold on %#3e7749');%'#3e7749'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,7000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(SS_neu_booster05_6m(:,10000)), median(SS_neu_booster05_6m(:,10000))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#3e7749');
hold on
plot([median(SS_neu_booster1(:,7000)), median(SS_neu_booster1(:,7000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');

xlim([20 100])
yticks(0:0.2:0.6)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_5_1y6m');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters');
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

figure(84)
h1 = histogram(SS_neu_booster05_6m(:,12000),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.3);
hold on %#3e7749');%'#3e7749'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,8000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(SS_neu_booster05_6m(:,12000)), median(SS_neu_booster05_6m(:,12000))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#3e7749');
hold on
plot([median(SS_neu_booster1(:,8000)), median(SS_neu_booster1(:,8000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
%ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.6)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_6_1y6m');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%
figure(85)
h1 = histogram(SS_neu_booster05_6m(:,14000),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.3);
hold on %#3e7749');%'#3e7749'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,9000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(SS_neu_booster05_6m(:,14000)), median(SS_neu_booster05_6m(:,14000))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#3e7749');
hold on
plot([median(SS_neu_booster1(:,9000)), median(SS_neu_booster1(:,9000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
%ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.6)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_7_1y6m');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

Edges = 0:10:100;
figure(86)
h1 = histogram(SS_neu_booster05_6m(:,16000),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.3);
hold on %#3e7749');%'#3e7749'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(SS_neu_booster1(:,10000),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(SS_neu_booster05_6m(:,16000)), median(SS_neu_booster05_6m(:,16000))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#3e7749');
hold on
plot([median(SS_neu_booster1(:,10000)), median(SS_neu_booster1(:,10000))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
%ylim([0 0.4])
xlim([20 100])
yticks(0:0.2:0.6)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'Neutralization_8_1y6m');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.2 2.5]); 
set(gcf, 'PaperSize', [3.2 2.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

close all
%
%%
m1 =[median(HH_neu_booster05(:,4000)),median(HH_neu_booster05(:,5000)),median(HH_neu_booster05(:,6000)),...
    median(HH_neu_booster05(:,7000)),median(HH_neu_booster05(:,8000)),median(HH_neu_booster05(:,9000)),...
    median(HH_neu_booster05(:,10000))];
m2 =[median(SS_neu_booster05(:,4000)),median(SS_neu_booster05(:,5000)),median(SS_neu_booster05(:,6000)),...
    median(SS_neu_booster05(:,7000)),median(SS_neu_booster05(:,8000)),median(SS_neu_booster05(:,9000)),...
    median(SS_neu_booster05(:,10000))];

m3 =[median(HH_neu_booster1(:,4000)),median(HH_neu_booster1(:,5000)),median(HH_neu_booster1(:,6000)),...
    median(HH_neu_booster1(:,7000)),median(HH_neu_booster1(:,8000)),median(HH_neu_booster1(:,9000)),...
    median(HH_neu_booster1(:,10000))];
m4 =[median(SS_neu_booster1(:,4000)),median(SS_neu_booster1(:,5000)),median(SS_neu_booster1(:,6000)),...
    median(SS_neu_booster1(:,7000)),median(SS_neu_booster1(:,8000)),median(SS_neu_booster1(:,9000)),...
    median(SS_neu_booster1(:,10000))];
m5 = [median(SS_neu_booster05_6m(:,4000)),median(SS_neu_booster05_6m(:,6000)),median(SS_neu_booster05_6m(:,8000)),...
    median(SS_neu_booster05_6m(:,10000)),median(SS_neu_booster05_6m(:,12000)),median(SS_neu_booster05_6m(:,14000)),...
    median(SS_neu_booster05_6m(:,16000))];

m = round([m1;m2;m3;m4;m5],2);
%%
b10 = sort(SS_neu_booster05(:,4000), 'descend'); 
b11 = sort(SS_neu_booster05(:,5000), 'descend'); 
b12 = sort(SS_neu_booster05(:,6000), 'descend'); 
b13 = sort(SS_neu_booster05(:,7000), 'descend'); 
b14 = sort(SS_neu_booster05(:,8000), 'descend'); 
b15 = sort(SS_neu_booster05(:,9000), 'descend'); 
b16 = sort(SS_neu_booster05(:,10000), 'descend'); 
b1 = [b10, b11,b12,b13,b14,b15,b16];

b20 = sort(SS_neu_booster1(:,4000), 'descend');
b21 = sort(SS_neu_booster1(:,5000), 'descend'); 
b22 = sort(SS_neu_booster1(:,6000), 'descend'); 
b23 = sort(SS_neu_booster1(:,7000), 'descend'); 
b24 = sort(SS_neu_booster1(:,8000), 'descend'); 
b25 = sort(SS_neu_booster1(:,9000), 'descend'); 
b26 = sort(SS_neu_booster1(:,10000), 'descend'); 
b2 = [b20, b21,b22,b23,b24,b25,b26];

b30 = sort(SS_neu_booster05_6m(:,4000), 'descend'); 
b31 = sort(SS_neu_booster05_6m(:,6000), 'descend'); 
b32 = sort(SS_neu_booster05_6m(:,8000), 'descend'); 
b33 = sort(SS_neu_booster05_6m(:,10000), 'descend'); 
b34 = sort(SS_neu_booster05_6m(:,12000), 'descend'); 
b35 = sort(SS_neu_booster05_6m(:,14000), 'descend'); 
b36 = sort(SS_neu_booster05_6m(:,16000), 'descend'); 
b3 = [b30, b31,b32,b33,b34,b35,b36];

a10 = sort(HH_neu_booster05(:,4000), 'descend');
a11 = sort(HH_neu_booster05(:,5000), 'descend'); 
a12 = sort(HH_neu_booster05(:,6000), 'descend'); 
a13 = sort(HH_neu_booster05(:,7000), 'descend'); 
a14 = sort(HH_neu_booster05(:,8000), 'descend'); 
a15 = sort(HH_neu_booster05(:,9000), 'descend'); 
a16 = sort(HH_neu_booster05(:,10000), 'descend'); 
a1 = [a10, a11,a12,a13,a14,a15,a16];

a20 = sort(HH_neu_booster1(:,4000), 'descend');
a21 = sort(HH_neu_booster1(:,5000), 'descend'); 
a22 = sort(HH_neu_booster1(:,6000), 'descend'); 
a23 = sort(HH_neu_booster1(:,7000), 'descend'); 
a24 = sort(HH_neu_booster1(:,8000), 'descend'); 
a25 = sort(HH_neu_booster1(:,9000), 'descend'); 
a26 = sort(HH_neu_booster1(:,10000), 'descend'); 
a2 = [a20, a21,a22,a23,a24,a25,a26];

