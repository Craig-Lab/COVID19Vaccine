clear all
format long
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
senior_LLL= load('neu_senior_LLL.mat').LLL_neu;
senior_LHL= load('neu_senior_LHL.mat').LHL_neu;
senior_LLH= load('neu_senior_LLH.mat').LLH_neu;
senior_LHH= load('neu_senior_LHH.mat').LHH_neu;
senior_HLL= load('neu_senior_HLL.mat').HLL_neu;
senior_HHL= load('neu_senior_HHL.mat').HHL_neu;
senior_HLH= load('neu_senior_HLH.mat').HLH_neu;
senior_HHH= load('neu_senior_HHH.mat').HHH_neu;

HCW_LLL= load('neu_HCW_LLL.mat').LLL_neu;
HCW_LHL= load('neu_HCW_LHL.mat').LHL_neu;
HCW_LLH= load('neu_HCW_LLH.mat').LLH_neu;
HCW_LHH= load('neu_HCW_LHH.mat').LHH_neu;
HCW_HLL= load('neu_HCW_HLL.mat').HLL_neu;
HCW_HHL= load('neu_HCW_HHL.mat').HHL_neu;
HCW_HLH= load('neu_HCW_HLH.mat').HLH_neu;
HCW_HHH= load('neu_HCW_HHH.mat').HHH_neu;

time = load('neu_HCW_HHH.mat').time;
s_LLL = median(senior_LLL);
s_LHL = median(senior_LHL);
s_LLH = median(senior_LLH);
s_LHH = median(senior_LHH);
s_HLL = median(senior_HLL);
s_HHL = median(senior_HHL);
s_HLH = median(senior_HLH);
s_HHH = median(senior_HHH);
h_LLL = median(HCW_LLL);
h_LHL = median(HCW_LHL);
h_LLH = median(HCW_LLH);
h_LHH = median(HCW_LHH);
h_HLL = median(HCW_HLL);
h_HHL = median(HCW_HHL);
h_HLH = median(HCW_HLH);
h_HHH = median(HCW_HHH);

%%
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');

map_senior ={ '#8CA6BD', '#6E95B7', '#465F75','#556573'};
%map_senior ={ '#495971', '#526683', '#5b7395', '#6580a8', '#6e8cba','#7799cc'};
col_days_senior = validatecolor(map_senior, 'multiple');
map_hcw ={'#d5abab',  '#b76767',  '#ac1c1e', '#8b0000'};

%map_hcw ={'#d5abab', '#c68989', '#b76767', '#a94545', '#9a2222', '#8b0000'};
col_days_hcw =validatecolor(map_hcw, 'multiple');
%map_hcw={'#fda057', '#e05206','#c87e7b','#8b0000'};
%map ={'#e9e9e9', '#d9d9d9', '#c6c6c6', '#b0b0b0', '#959595','#7e7e7e', '#686868', '#515151', '#333333', '#181818'};
map= {'#F8F9FA','#E9ECEF','#DEE2E6','#CED4DA','#ADB5BD'};
col_days = validatecolor(map, 'multiple');
x = [0 700 700 0];
y1 =[70 70 60 60] ;
y2 = [80 80 70 70];
y3 =  [90 90 80 80];
y4 = [95 95 90 90];
y5 = [100 100 95 95];
figure(1)
patch(x,y1,col_days(1,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y2,col_days(2,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y3,col_days(3,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y4,col_days(4,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y5,col_days(5,:),'LineStyle','none','FaceAlpha', 1);
hold on
plot (time,s_HLL,'--','Color',col_days_senior(1,:),'LineWidth',1.7)
hold on
plot (time,s_HLH,'--','Color',col_days_senior(2,:),'LineWidth',1.7)
hold on
plot (time,s_HHL,'--','Color',col_days_senior(3,:),'LineWidth',1.7)
hold on
plot (time,s_HHH,'-','Color',col_days_senior(4,:),'LineWidth',1.5)
hold on
plot (time,h_HLL,'--','Color',col_days_hcw(1,:),'LineWidth',1.7)
hold on
plot (time,h_HLH,'--','Color',col_days_hcw(2,:),'LineWidth',1.7)
hold on
plot (time,h_HHL,'-','Color',col_days_hcw(3,:),'LineWidth',1.5)
hold on
plot (time,h_HHH,'--','Color',col_days_hcw(4,:),'LineWidth',1.7)
hold on
xlim([0 630])
xticks([0,60,270,630])
box off
ylim([0 100])
xlim([0 630])
ylabel( 'Neutralization (%)')

xlabel('Time (Days) ')
figfile = fullfile(pathname,'DoseSize_median');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 8 4.5]); 
set(gcf, 'PaperSize', [8 4.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%%
figure(3)
patch(x,y1,col_days(1,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y2,col_days(2,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y3,col_days(3,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y4,col_days(4,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y5,col_days(5,:),'LineStyle','none','FaceAlpha', 1);
hold on
h1 = plot(time, s_LLL, 'Color', col_days_senior(1,:), 'LineWidth', 1.5);
hold on
h2 = plot(time, s_LLH, 'Color',col_days_senior(2,:), 'LineWidth', 1.5);
hold on
h3 = plot(time, s_LHL, 'Color', col_days_senior(3,:), 'LineWidth', 1.5);
hold on
h4 = plot(time, s_LHH, 'Color',col_days_senior(4,:), 'LineWidth', 1.5);
hold on
h5 = plot(time, s_HLL, '--', 'Color',col_days_hcw(1,:), 'LineWidth', 1.5);
hold on
h6 = plot(time, s_HLH, '--', 'Color',col_days_hcw(2,:), 'LineWidth', 1.5);
hold on
h7 = plot(time, s_HHL, '--', 'Color', col_days_hcw(3,:), 'LineWidth', 1.5);
hold on
h8 = plot(time, s_HHH, '--', 'Color', col_days_hcw(4,:), 'LineWidth', 1.5);
hold off


xlim([0 630])
xticks([0, 60, 270, 630])
box off
ylim([0 100])
xlim([0 630])
ylabel('Neutralization (%)')

xlabel('Time (Days) ')
figfile = fullfile(pathname,'DoseSize_median_seniors');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gcf, 'PaperSize', [7.5 4.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%
figure(2)
patch(x,y1,col_days(1,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y2,col_days(2,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y3,col_days(3,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y4,col_days(4,:),'LineStyle','none','FaceAlpha', 1);
hold on
patch(x,y5,col_days(5,:),'LineStyle','none','FaceAlpha', 1);
hold on
h1 = plot(time, h_LLL, 'Color', col_days_senior(1,:), 'LineWidth', 1.5);
hold on
h2 = plot(time, h_LLH, 'Color',col_days_senior(2,:), 'LineWidth', 1.5);
hold on
h3 = plot(time, h_LHL, 'Color', col_days_senior(3,:), 'LineWidth', 1.5);
hold on
h4 = plot(time, h_LHH, 'Color',col_days_senior(4,:), 'LineWidth', 1.5);
hold on
h5 = plot(time, h_HLL, '--', 'Color',col_days_hcw(1,:), 'LineWidth', 1.5);
hold on
h6 = plot(time, h_HLH, '--', 'Color',col_days_hcw(2,:), 'LineWidth', 1.5);
hold on
h7 = plot(time, h_HHL, '--', 'Color', col_days_hcw(3,:), 'LineWidth', 1.5);
hold on
h8 = plot(time, h_HHH, '--', 'Color', col_days_hcw(4,:), 'LineWidth', 1.5);
hold off

xlim([0 630])
xticks([0, 60, 270, 630])
box off
ylim([0 100])
xlim([0 630])
ylabel('Neutralization (%)')

xlabel('Time (Days) ')
figfile = fullfile(pathname,'DoseSize_median_hcws');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gcf, 'PaperSize', [7.5 4.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%
h3_hcw = h_HHH(2001:4000);
subset_h3_hcw = h3_hcw(250:2000);
index_in_subset = find(subset_h3_hcw < 90, 1);
hcw_90= time(2250 + index_in_subset - 1)-270;

h3_senior = s_HHH(2001:4000);
subset_h3_senior = h3_senior(250:2000);
index_in_subset = find(subset_h3_senior < 90, 1);
senior_90 = time(2250 + index_in_subset - 1)-270;
hcw_90-senior_90

median(h_HHH(4000))
median(s_HHH(4000))