clear all
format long
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
senior_int12_1m= load('Neu_senior_int12_1m.mat').SH;
senior_int12_2m= load('Neu_senior_int12_2m.mat').SH;
senior_int12_4m= load('Neu_senior_int12_4m.mat').SH;
senior_int12_6m= load('Neu_senior_int12_6m.mat').SH;
senior_int23_6m= load('Neu_senior_int23_6m.mat').SH;
senior_int23_7m= load('Neu_senior_int23_7m.mat').SH;
senior_int23_8m= load('Neu_senior_int23_8m.mat').SH;
senior_int23_9m= load('Neu_senior_int23_9m.mat').SH;
time1 = load('Neu_senior_int12_1m.mat').time;
time2= load('Neu_senior_int12_2m.mat').time;
time3 = load('Neu_senior_int12_4m.mat').time;
time4 = load('Neu_senior_int12_6m.mat').time;
time5 = load('Neu_senior_int23_6m.mat').time;
time6 = load('Neu_senior_int23_7m.mat').time;
time7 = load('Neu_senior_int23_8m.mat').time;
time8 = load('Neu_senior_int23_9m.mat').time;
hcw_int12_1m= load('Neu_hcw_int12_1m.mat').SH;
hcw_int12_2m= load('Neu_hcw_int12_2m.mat').SH;
hcw_int12_4m= load('Neu_hcw_int12_4m.mat').SH;
hcw_int12_6m= load('Neu_hcw_int12_6m.mat').SH;
hcw_int23_6m= load('Neu_hcw_int23_6m.mat').SH;
hcw_int23_7m= load('Neu_hcw_int23_7m.mat').SH;
hcw_int23_8m= load('Neu_hcw_int23_8m.mat').SH;
hcw_int23_9m= load('Neu_hcw_int23_9m.mat').SH;

median_senior_int12_1m = median(senior_int12_1m, 1);
median_senior_int12_2m = median(senior_int12_2m, 1);
median_senior_int12_4m = median(senior_int12_4m, 1);
median_senior_int12_6m = median(senior_int12_6m, 1);
median_senior_int23_6m = median(senior_int23_6m, 1);
median_senior_int23_7m = median(senior_int23_7m, 1);
median_senior_int23_8m = median(senior_int23_8m, 1);
median_senior_int23_9m = median(senior_int23_9m, 1);

median_hcw_int12_1m = median(hcw_int12_1m, 1);
median_hcw_int12_2m = median(hcw_int12_2m, 1);
median_hcw_int12_4m = median(hcw_int12_4m, 1);
median_hcw_int12_6m = median(hcw_int12_6m, 1);
median_hcw_int23_6m = median(hcw_int23_6m, 1);
median_hcw_int23_7m = median(hcw_int23_7m, 1);
median_hcw_int23_8m = median(hcw_int23_8m, 1);
median_hcw_int23_9m = median(hcw_int23_9m, 1);
%%
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');

map_senior ={ '#8CA6BD', '#5B7B98', '#465F75','#556573'};
%map_senior ={ '#495971', '#526683', '#5b7395', '#6580a8', '#6e8cba','#7799cc'};
col_days_senior = validatecolor(map_senior, 'multiple');
map_hcw ={'#d5abab',  '#b76767',  '#ac1c1e', '#8b0000'};

%map_hcw ={'#d5abab', '#c68989', '#b76767', '#a94545', '#9a2222', '#8b0000'};
col_days_hcw =validatecolor(map_hcw, 'multiple');
map_hcw={'#fda057', '#e05206','#c87e7b','#8b0000'};
%map ={'#e9e9e9', '#d9d9d9', '#c6c6c6', '#b0b0b0', '#959595','#7e7e7e', '#686868', '#515151', '#333333', '#181818'};
map= {'#F8F9FA','#E9ECEF','#DEE2E6','#CED4DA','#ADB5BD'};
col_days = validatecolor(map, 'multiple');
x = [0 700 700 0];
y1 =[70 70 60 60] ;
y2 = [80 80 70 70];
y3 =  [90 90 80 80];
y4 = [95 95 90 90];
y5 = [100 100 95 95];
%%
figure(1)
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
plot (time1,median_senior_int12_1m,'Color',col_days_senior(1,:),'LineWidth',1.2)
hold on
plot (time2,median_senior_int12_2m,'-.','Color',col_days_senior(2,:),'LineWidth',1.2)
hold on
plot (time3,median_senior_int12_4m,'-','Color',col_days_senior(3,:),'LineWidth',1.2)
hold on
plot (time4,median_senior_int12_6m,'-','Color',col_days_senior(4,:),'LineWidth',1.2)
hold off
xlim([0 630])
xticks([0,60,270,630])
box off
ylim([20 100])
xlim([0 630])
ylabel( 'Neutralization (%)')

xlabel('Time (Days) ')
figfile = fullfile(pathname,'Senior_int12');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gcf, 'PaperSize', [7.5 4.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
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
plot (time5,median_senior_int23_6m,'Color',col_days_senior(1,:),'LineWidth',1.2)
hold on
plot (time6,median_senior_int23_7m,'-.','Color',col_days_senior(2,:),'LineWidth',1.2)
hold on
plot (time7,median_senior_int23_8m,'-','Color',col_days_senior(3,:),'LineWidth',1.2)
hold on
plot (time8,median_senior_int23_9m,'-','Color',col_days_senior(4,:),'LineWidth',1.2)
hold off
xlim([0 630])
xticks([0,60,270,630])
box off
ylim([20 100])
xlim([0 630])
ylabel( 'Neutralization (%)')

xlabel('Time (Days) ')
figfile = fullfile(pathname,'Senior_int23');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gcf, 'PaperSize', [7.5 4.5]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
median_senior_int23_9m(end)
median_senior_int12_4m(end)
