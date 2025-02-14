clear all

clear all
senior_par_data = readtable('senior_ant_par.xlsx');
senior_dt= table2array(senior_par_data(:,2));
senior_lam2 = table2array(senior_par_data(:,3));
senior_lam3 = table2array(senior_par_data(:,4));
senior_age = table2array(senior_par_data(:,5));

hcw_par_data = readtable('hcw_ant_par.xlsx');
hcw_dt= table2array(hcw_par_data(:,2));
hcw_lam2 = table2array(hcw_par_data(:,3));
hcw_lam3 = table2array(hcw_par_data(:,4));
hcw_age = table2array(hcw_par_data(:,5));

hcw_neu_data = readtable('hcw_neu_par.xlsx');
hcw_E0= table2array(hcw_neu_data(:,2));
hcw_EC50 = table2array(hcw_neu_data(:,3));
hcw_h = table2array(hcw_neu_data(:,4));

senior_neu_data = readtable('senior_neu_par.xlsx');
senior_E0= table2array(senior_neu_data(:,2));
senior_EC50 = table2array(senior_neu_data(:,3));
senior_h = table2array(senior_neu_data(:,4));

H = [hcw_age, hcw_dt, hcw_lam2, hcw_lam3, hcw_E0,hcw_EC50, hcw_h];%, hcw_par];
S= [senior_age, senior_dt, senior_lam2, senior_lam3, senior_E0,senior_EC50, senior_h];%, sen_par];
A = [H;S];
AA = sortrows(A,5);
[R,P] = corrcoef(A,'Alpha',0.05);
a = round([min(H);max(H);median(H);std(H)],2);%HCW
b = round([min(S);max(S);median(S); std(S)],2);%Senior
a_median= median(H);
b_median =median(S);
aa = [a_median;b_median];
%%

lam_hcw = [hcw_lam2, hcw_lam3]; % Array for HCW
lam_senior = [senior_lam2, senior_lam3]; % Array for Seniors

% Calculate cumulative sums
lam_sum_hcw =[1+median(hcw_lam2), 1+median(hcw_lam2+hcw_lam3)];
lam_sum_senior =[1+median(senior_lam2), 1+median(senior_lam2+senior_lam3)];
total_hcw = 1+hcw_lam2+hcw_lam3;
total_senior = 1+senior_lam2+senior_lam3;
aaa= [min(total_hcw),min(total_senior);max(total_senior),max(total_senior);...
    std(total_hcw),std(total_senior);median(total_hcw),median(total_senior)];


%%
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');

% Define time intervals
t = [0, 60, 60, 270, 270, 630, 630, 640]; % Time points
value1 = 1;                               % Constant value for t < 60

% Healthcare workers' data
value21_hcw = 1 + hcw_lam2;                          % Amplification at t=270
value31_hcw = 1 + hcw_lam2 + hcw_lam3;               % Amplification at t=630

% Seniors' data
value22_senior = 1 + senior_lam2;                    % Amplification at t=270
value32_senior = 1 + senior_lam2 + senior_lam3;      % Amplification at t=630

%Kolmogorov-Smirnov test
[h1, p1, ksstat1] = kstest2(value21_hcw, value22_senior);
[h2, p2, ksstat2] = kstest2(value31_hcw, value32_senior);
ks = [h1,h2;p1,p2];

% Calculate medians
median_value21_hcw = median(value21_hcw);
median_value31_hcw = median(value31_hcw);
median_value22_senior = median(value22_senior);
median_value32_senior = median(value32_senior);

% Calculate 95% confidence intervals
CI95_hcw2 = prctile(value21_hcw, [2.5, 97.5]);
CI95_hcw3 = prctile(value31_hcw, [2.5, 97.5]);
CI95_senior2 = prctile(value22_senior, [2.5, 97.5]);
CI95_senior3 = prctile(value32_senior, [2.5, 97.5]);

%%  Plot for HCWs
figure;
%95% CI for HCWs
x_CI_hcw = [60, 270, 270, 630, 630, 270, 270, 60];
y_CI_hcw = [CI95_hcw2(1), CI95_hcw2(1), CI95_hcw3(1), CI95_hcw3(1), CI95_hcw3(2), CI95_hcw3(2), CI95_hcw2(2), CI95_hcw2(2)];
patch(x_CI_hcw, y_CI_hcw, 1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.2);
hold on;

% ---- 95% CI for Seniors ----
x_CI_senior = [60, 270, 270, 630, 630, 270, 270, 60];
y_CI_senior = [CI95_senior2(1), CI95_senior2(1), CI95_senior3(1), CI95_senior3(1), CI95_senior3(2), CI95_senior3(2), CI95_senior2(2), CI95_senior2(2)];
patch(x_CI_senior, y_CI_senior,1, 'facecolor','#4276bc','edgecolor','none','facealpha', 0.2); %ci;

% ---- Median Step Function for HCWs ----
values_median_hcw = [value1, value1, median_value21_hcw, median_value21_hcw, median_value31_hcw, median_value31_hcw, median_value31_hcw, median_value31_hcw];
stairs(t, values_median_hcw, 'Color', '#b24a4b', 'LineWidth', 1.2, 'DisplayName', 'Median HCW');

% ---- Median Step Function for Seniors ----
values_median_senior = [value1, value1, median_value22_senior, median_value22_senior, median_value32_senior, median_value32_senior, median_value32_senior, median_value32_senior];
stairs(t, values_median_senior,  'LineWidth', 1.2, 'Color', '#4276bc','DisplayName', 'Median Senior');


% Customize plot
ylabel('Immunol amplification');
xlabel('Time');
xlim([0, 630]);
ylim([0, max([CI95_hcw3(2), CI95_senior3(2)]) + 1]);
xticks([0, 60, 270, 600]);
ylim([0 30])
yticks(0:10:30)
figfile = fullfile(pathname,'ImmAmp');

set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',7);

set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%
figure(2)
h1 = histogram(value21_hcw,0:2:30, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.2);
hold on 
h2 = histogram(value22_senior,0:2:30, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.2);
hold on
plot([median(value21_hcw), median(value21_hcw)],[0 h1.Values(5)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(value22_senior), median(value22_senior)],[0 h2.Values(5)],'LineWidth',1.5, 'Color','#4276bc');
ylim([0 0.6])
xlim([4 16])
yticks(0:0.2:0.6)
ylabel('Probability')
xlabel('Immunological amplification');
box off
figfile = fullfile(pathname,'ImmunAmp2');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',7);

set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%%
figure(3)
h1 = histogram(value31_hcw,0:2:30, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.2);
hold on 
h2 = histogram(value32_senior,0:2:30, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.2);
hold on
plot([median(value31_hcw), median(value31_hcw)],[0 h1.Values(9)],'LineWidth',1.2, 'Color','#b24a4b');
hold on
plot([median(value32_senior), median(value32_senior)],[0 h2.Values(9)],'LineWidth',1.2, 'Color','#4276bc');
ylim([0 0.3])
xlim([10 28])
yticks(0:0.1:0.3)
ylabel('Probability')
xlabel('Immunological amplification');
box off
figfile = fullfile(pathname,'ImmunAmpTotal');

set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',7);

set(gcf, 'PaperPosition', [0 0 6.5 3.3]); 
set(gcf, 'PaperSize', [6.5 3.3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);


