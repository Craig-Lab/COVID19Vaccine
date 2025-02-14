clear all

Anti_data = readtable('all.csv');
Neut_data = readtable('neutralization.csv');

Age = categorical(Anti_data.Type);
idx_senior = (Age=='SENIOR');
idx_hcw = (Age == 'HCW');

senior_data=table2array(Neut_data(idx_senior,2:7));
hcw_data=table2array(Neut_data(idx_hcw,2:7));
median_NS= median(senior_data,"omitnan");
std_NS = std(senior_data,"omitnan");
median_NH= median(hcw_data,"omitnan");
std_NH =  std(hcw_data,"omitnan");

hcwdata = readtable('hcw.csv');
seniordata = readtable('senior.csv');
Anti_H = table2array(hcwdata(:,2:7));
Anti_S = table2array(seniordata(:,2:7));

AH = Anti_H;
AS =  Anti_S;

median_AS =median(AS/1000,"omitnan");%predicted antibody level of senior
median_AH = median(AH/1000,"omitnan"); %predicted antibody level of hcw

[rho1,pval1] = corr(median_NS', median_AS');
[rho2,pval2] = corr(median_NH', median_AH');
%%
Timepoints= [30, 90, 150,240,300,360];
senior_ci_data_raw = load('senior_neu_ind.mat');
senior_ci_data = senior_ci_data_raw.SS;
SS_neu_ci = zeros(2,4000);
hcw_ci_data_raw = load('hcw_neu_ind.mat');
hcw_ci_data = hcw_ci_data_raw.HH;
HH_neu_ci = zeros(2,4000);
N1 = median(hcw_ci_data);
N2 = median(senior_ci_data);
time = senior_ci_data_raw.Time;
for i = 1:4000
    hcw_data = hcw_ci_data(:,i);
    senior_data = senior_ci_data(:,i);    

    HH_neu_ci(1,i) = prctile(hcw_data , 97.5); % Upper bound (97.5th percentile)
    HH_neu_ci(2,i) = prctile(hcw_data , 2.5); % Lower bound (2.5th percentile)
    
    SS_neu_ci(1,i) = prctile(senior_data, 97.5); % Upper bound (97.5th percentile)
    SS_neu_ci(2,i) = prctile(senior_data, 2.5); % Lower bound (2.5th percentile)
end
map= {'#F8F9FA','#E9ECEF','#DEE2E6','#CED4DA','#ADB5BD'};
col_days = validatecolor(map, 'multiple');
x = [0 700 700 0];
y1 =[70 70 60 60] ;
y2 = [80 80 70 70];
y3 =  [90 90 80 80];
y4 = [95 95 90 90];
y5 = [100 100 95 95];
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');

figure(4)
% patch(x,y1,col_days(1,:),'LineStyle','none')
% hold on
% patch(x,y2,col_days(2,:),'LineStyle','none')
% hold on
% patch(x,y3,col_days(3,:),'LineStyle','none')
% hold on
% patch(x,y4,col_days(4,:),'LineStyle','none')
% hold on
% patch(x,y5,col_days(5,:),'LineStyle','none')
idxDashed = time >= 0 & time <= 60;
idxSolid = time > 60;
patch([time, fliplr(time)],[SS_neu_ci(1,:), fliplr(SS_neu_ci(2,:))],1,'facecolor','#4276bc','edgecolor','none','facealpha', 0.15); %ci
hold on
patch([time, fliplr(time)],[HH_neu_ci(1,:), fliplr(HH_neu_ci(2,:))],1,'facecolor','#b24a4b','edgecolor','none','facealpha', 0.15); %ci
hold on
plot(time(idxDashed), N1(idxDashed), 'LineWidth', 1.2, 'Color', '#b24a4b', 'LineStyle', '-.');
hold on;
plot(time(idxSolid), N1(idxSolid), 'LineWidth', 1.2, 'Color', '#b24a4b');
hold on;
plot(time(idxDashed), N2(idxDashed), 'LineWidth', 1.2, 'Color', '#4276bc', 'LineStyle', '-.');
hold on;
plot(time(idxSolid), N2(idxSolid), 'LineWidth', 1.2, 'Color', '#4276bc');
hold on;
%errorbar(Timepoints, median_NS,std_NS,'o','MarkerEdgeColor','#4276bc','MarkerFaceColor','#4276bc','MarkerSize',5,'Capsize',5,'LineWidth',1.5,'Color','#4276bc');
hold on
%errorbar(Timepoints, median_NH,std_NH,'o','MarkerEdgeColor','#b24a4b','MarkerFaceColor','#b24a4b','MarkerSize',5,'Capsize',5,'LineWidth',1.5,'Color','#b24a4b');
hold on
h3 = scatter(Timepoints, median_NH,15, 'o','filled','MarkerFaceColor','#b24a4b');
hold on
h4 = scatter(Timepoints,median_NS,15, 'o','filled','MarkerFaceColor','#4276bc');
hold off
xlim([0 630])
ylim([0 100])
xlabel('Time (days)')
ylabel('Neutralization (%)')

box off
figfile = fullfile(pathname,'Predicted_Neutralization');

set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',8);
set(gcf, 'PaperPosition', [0 0 6.8 4.2]); 
set(gcf, 'PaperSize', [6.8 4.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);


%% Neutralization duriblity
threshold = [60, 70, 80, 90, 95];
I1= zeros(32,5); % hcw_ci_data
I2= zeros(27,5); % senior_ci_data
I3= zeros(32,5); % senior_ci_data
I4= zeros(27,5); % senior_ci_data
[a1, IP1] = max(hcw_ci_data(:,1001:2000),[],2); %hcw_ci_data
[a2, IP2] = max(senior_ci_data(:,1001:2000),[],2); %senior_ci_data
[a3, IP3] = max(hcw_ci_data(:,2001:4000),[],2); %hcw_ci_data
[a4, IP4] = max(senior_ci_data(:,2001:4000),[],2); %senior_ci_data

IP1 = IP1+1000;
IP2 = IP2+1000;
IP3 = IP3+2000;
IP4 = IP4+2000;

for i = 1:32
    for j = 1:5
            idx = IP1(i);
                if hcw_ci_data(i,idx)<threshold(j)
                    I1(i,j) =1001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(hcw_ci_data(i,idx:2000)>threshold(j)));
                     I1(i,j)= a+idx-1;
                end
    end
end
for i = 1:27
    for j = 1:5
            idx = IP2(i);
                if senior_ci_data(i,idx)<threshold(j)
                    I2(i,j) =1001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(senior_ci_data(i,idx:2000)>threshold(j)));
                     I2(i,j)= a+idx-1;
                end
    end
end
for i = 1:32
    for j = 1:5
            idx = IP3(i);
                if hcw_ci_data(i,idx)<threshold(j)
                    I3(i,j) =2001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(hcw_ci_data(i,idx:4000)>threshold(j)));
                     I3(i,j)= a+idx-1;
                end
    end
end

for i = 1:27
    for j = 1:5
            idx = IP4(i);
                if senior_ci_data(i,idx)<threshold(j)
                    I4(i,j) =2001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(senior_ci_data(i,idx:4000)>threshold(j)));
                     I4(i,j)= a+idx-1;
                end
    end
end
T1 = round(time(I1))-60;
T2 = round(time(I2))-60;
T3 = round(time(I3));
T4 = round(time(I4));
T3(T3 >= 629)= 630;
T4(T4 >= 629)= 630;
T4 = T4-270;
T3 = T3-270;

m1 =median(T1);
m2 = median(T2);
m3 = median(T3);
m4 = median(T4);

m = [m1;m2;m3;m4];
min_a = [min(a1),min(a2),min(a3),min(a4)];
max_a = [max(a1),max(a2),max(a3),max(a4)];
median_a = [median(a1),median(a2),median(a3),median(a4)];
a_ana = [min_a;max_a;median_a];
min_t= [min(T1);min(T2);min(T3);min(T4)]';
max_t = [max(T1);max(T2);max(T3);max(T4)]';
interval_ana = [m1',m2',m3',m4'];
%% analyzing max neutralization
[a11,b11] = max(N1(1:2000));
peak1_hcw_sd = std(hcw_ci_data(:,b11));
[a12,b12] = max(N2(1:2000));
peak1_senior_sd = std(senior_ci_data(:,b12));
N_2dose = [hcw_ci_data(:,b11);senior_ci_data(:,b12)];
[a21,b21] = max(N1);
peak2_hcw_sd = std(hcw_ci_data(:,b21));
[a22,b22] = max(N2);
peak2_senior_sd = std(senior_ci_data(:,b22));
count = sum(N_2dose < 95);
T80 = sort([T3(:,3);T4(:,3)],'descend');
T90 = sort([T3(:,4);T4(:,4)],'descend');
T95 = sort([T3(:,5);T4(:,5)],'descend');

sum(T3(:,3)>200)%HCW >200 days (80% threshold)
sum(T4(:,3)>200)%Senior >200 days (80% threshold)
sum(T4(:,1)==360)
%%
senior_compare = [senior_ci_data(:,1000),senior_ci_data(:,1287),senior_ci_data(:,2000),senior_ci_data(:,2335),senior_ci_data(:,3168),senior_ci_data(:,4000)];
hcw_compare =[hcw_ci_data(:,1000),hcw_ci_data(:,1287),hcw_ci_data(:,2000),hcw_ci_data(:,2335),hcw_ci_data(:,3168),hcw_ci_data(:,4000)];
m10 = median(senior_compare);
m11 = median(hcw_compare);
mm =[m11;m10];


aa1= sort(senior_ci_data(:,4000),'descend');
aa2= sort(hcw_ci_data(:,4000),'descend');
aa3= sort(hcw_ci_data(:,2000),'descend');
hh =[aa3,aa2];
%% analyzing T matrix
%HCW
min_T1 = min(T1);
max_T1 = max(T1);
median_T1 = median(T1);
%senior
min_T2 = min(T2);
max_T2 = max(T2);
median_T2 = median(T2);
%HCW
min_T3 = min(T3);
max_T3 = max(T3);
median_T3 = median(T3);
%Senior
min_T4 = min(T4);
max_T4 = max(T4);
median_T4 = median(T4);

% Construct a summary matrix
SummaryMatrix = [
    min_T1, max_T1, median_T1;
    min_T2, max_T2, median_T2;
    min_T3, max_T3, median_T3;
    min_T4, max_T4, median_T4
];
%%
Edges = 0:30:390;
figure(31)
h1 = histogram(T1(:,1),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on 
h2 = histogram(T2(:,1),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
% Calculate cumulative probabilities and locate 50% for T1
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');
xlim([0 380])
ylim([0 0.4])
%xlim([0 220])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time2to3_60_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

%%

figure(32)
h1 = histogram(T1(:,2),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on 
h2 = histogram(T2(:,2),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');
ylim([0 0.4])
%xlim([0 220])
xlim([0 380])
yticks(0:0.1:0.4)
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time2to3_70_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI close(gcf);
close(gcf);



%%
figure(33)
h1 = histogram(T1(:,3),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3)
hold on
h2 = histogram(T2(:,3),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3)
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
hold on
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');
ylim([0 0.4])
xlim([0 380])
%xlim([0 220])
yticks(0:0.2:0.4)
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time2to3_80_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);


%%
figure(34)
h1 = histogram(T1(:,4),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(T2(:,4),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
hold on
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');
ylim([0 0.6])
%xlim([0 220])
%yticks(0:0.1:0.5)
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time2to3_90_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);


%%
figure(35)
h1 = histogram(T1(:,5),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(T2(:,5),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
hold on
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');
%ylim([0 0.7])
%xlim([0 220])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Time2to3_95_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

%%
Edges = 0:30:390;
figure(41)
h1 = histogram(T3(:,1),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on 
h2 = histogram(T4(:,1),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
hold on
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');
ylim([0 0.7])
yticks([0, 0.2,0.4, 0.6,0.7])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_60_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);


%%
figure(42)
h1 = histogram(T3(:,2),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on 
h2 = histogram(T4(:,2),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
hold on
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');ylim([0 0.5])
xlim([0 380])
ylim([0 0.6])
yticks(0:0.2:0.6)
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_70_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);


%%
figure(43)
h1 = histogram(T3(:,3),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3)
hold on
h2 = histogram(T4(:,3),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3)
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
hold on
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');hold off
ylim([0 0.3])
yticks(0:0.1:0.3)
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_80_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

%%
figure(44)
h1 = histogram(T3(:,4),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(T4(:,4),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
hold on
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');hold on
ylim([0 0.3])
xlim([0 380])
yticks([0, 0.1, 0.2,0.3])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_90_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

%%
figure(45)
h1 = histogram(T3(:,5),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(T4(:,5),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
cumulativeProb1 = cumsum(h1.Values);
binEdges1 = h1.BinEdges(1:end-1); % Left edges of bins
binIdx1 = find(cumulativeProb1 >= 0.5, 1); % Bin index where 50% is reached
binAt50_1 = binEdges1(binIdx1); % Bin edge at 50%
binHeight1 = h1.Values(binIdx1); % Height of the bin at 50%

% Calculate cumulative probabilities and locate 50% for T2
cumulativeProb2 = cumsum(h2.Values);
binEdges2 = h2.BinEdges(1:end-1); % Left edges of bins
binIdx2 = find(cumulativeProb2 >= 0.5, 1); % Bin index where 50% is reached
binAt50_2 = binEdges2(binIdx2); % Bin edge at 50%
binHeight2 = h2.Values(binIdx2); % Height of the bin at 50%

% Plot vertical lines at 50% cumulative probability with heights matching the bins
plot([binAt50_1, binAt50_1], [0, binHeight1], 'LineWidth', 1.5, 'Color', '#b24a4b');
hold on
plot([binAt50_2, binAt50_2], [0, binHeight2], 'LineWidth', 1.5, 'Color', '#4276bc');
ylim([0 0.3])
xlim([0 380])
yticks([0, 0.1, 0.2,0.3])
%ylabel('Probability')
xlabel('Time (days)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Time3to4_95_PoP');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

%%
Edges_new = 0:10:100;
figure(71)
h1 = histogram(hcw_ci_data(:,1000),Edges_new, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(senior_ci_data(:,1000),Edges_new, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(hcw_ci_data(:,1000)), median(hcw_ci_data(:,1000))],[0 h1.Values(4)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(senior_ci_data(:,1000)), median(senior_ci_data(:,1000))],[0 h2.Values(3)],'LineWidth',1.5,'Color','#4276bc');
hold on
ylim([0 0.8])
xlim([0 100])
%ylabel('Probability')
xlabel('Neutralization (%)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Neutralization_1st');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);




%%
[h1, p1, ksstat1] = kstest2(hcw_ci_data(:,1000), senior_ci_data(:,1000));
[h2, p2, ksstat2] = kstest2(hcw_ci_data(:,2000), senior_ci_data(:,2000));
[h3, p3, ksstat3] = kstest2(hcw_ci_data(:,4000), senior_ci_data(:,4000));
ks = [h1,h2,h3;p1,p2,p3];
[p4, h4, stats4] = ranksum(hcw_ci_data(:,1000), senior_ci_data(:,1000));
[p5, h5, stats5] = ranksum(hcw_ci_data(:,2000), senior_ci_data(:,2000));
[p6, h6, stats6] = ranksum(hcw_ci_data(:,4000), senior_ci_data(:,4000));
wmw = [h4,h5,h6;p4,p5,p6];


%%
Edges_new = 0:10:100;
figure(75)
h1 = histogram(hcw_ci_data(:,1287),Edges_new, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(senior_ci_data(:,1287),Edges_new, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(hcw_ci_data(:,1287)), median(hcw_ci_data(:,1287))],[0 h1.Values(9)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(senior_ci_data(:,1287)), median(senior_ci_data(:,1287))],[0 h2.Values(8)],'LineWidth',1.5,'Color','#4276bc');
hold on
ylim([0 0.4])
xlim([0 100])
%ylabel('Probability')
yticks([0 0.2 0.4])
xlabel('Neutralization (%)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Neutralization_120');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

%%
Edges_new = 0:10:100;
figure(75)
h1 = histogram(hcw_ci_data(:,2335),Edges_new, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(senior_ci_data(:,2335),Edges_new, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(hcw_ci_data(:,2335)), median(hcw_ci_data(:,2335))],[0 h1.Values(10)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(senior_ci_data(:,2335)), median(senior_ci_data(:,2335))],[0 h2.Values(10)],'LineWidth',1.5,'Color','#4276bc');
hold on
%ylim([0 0.5])
xlim([0 100])
%ylabel('Probability')
xlabel('Neutralization (%)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Neutralization_330');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);


%%
Edges_new = 0:10:100;
figure(75)
h1 = histogram(hcw_ci_data(:,1000),Edges_new, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(senior_ci_data(:,1000),Edges_new, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(hcw_ci_data(:,1000)), median(hcw_ci_data(:,1000))],[0 h1.Values(4)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(senior_ci_data(:,1000)), median(senior_ci_data(:,1000))],[0 h2.Values(3)],'LineWidth',1.5,'Color','#4276bc');
hold on
ylim([0 0.8])
xlim([0 100])
yticks(0:0.2:0.8)
%ylabel('Probability')
%xlabel('Neutralization (%)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Neutralization_60');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 2.5 3]); 
set(gcf, 'PaperSize', [2.5 3]);
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);


%
Edges_new = 0:10:100;
figure(73)
h1 = histogram(hcw_ci_data(:,2000),Edges_new, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(senior_ci_data(:,2000),Edges_new, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(hcw_ci_data(:,2000)), median(hcw_ci_data(:,2000))],[0 h1.Values(5)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(senior_ci_data(:,2000)), median(senior_ci_data(:,2000))],[0 h2.Values(4)],'LineWidth',1.5,'Color','#4276bc');
hold on
ylim([0 0.8])
xlim([0 100])
yticks(0:0.2:0.8)
%ylabel('Probability')
%xlabel('Neutralization (%)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Neutralization_2nd');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 2.5 3]); 
set(gcf, 'PaperSize', [2.5 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

Edges_new = 0:10:100;
figure(75)
h1 = histogram(hcw_ci_data(:,4000),Edges_new, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(senior_ci_data(:,4000),Edges_new, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([median(hcw_ci_data(:,4000)), median(hcw_ci_data(:,4000))],[0 h1.Values(7)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([median(senior_ci_data(:,4000)), median(senior_ci_data(:,4000))],[0 h2.Values(7)],'LineWidth',1.5,'Color','#4276bc');
hold on
ylim([0 0.8])
xlim([0 100])
yticks(0:0.2:0.8)
%ylabel('Probability')
%xlabel('Neutralization (%)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')
figfile = fullfile(pathname,'Neutralization_3rd');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 2.5 3]); 
set(gcf, 'PaperSize', [2.5 3]);
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf)


%%

Edges_new = 0:5:100;
figure(66)
h1 = histogram(max(hcw_ci_data(:,1001:2000), [], 2),Edges_new, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(max(senior_ci_data(:,1001:2000), [], 2),Edges_new, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([mean(max(hcw_ci_data(:,1001:2000), [], 2)), mean(max(hcw_ci_data(:,1001:2000), [], 2))],[0 h1.Values(19)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([mean(max(senior_ci_data(:,1001:2000), [], 2)), mean(max(senior_ci_data(:,1001:2000), [], 2))],[0 h2.Values(18)],'LineWidth',1.5,'Color','#4276bc');
hold on
ylim([0 0.6])
xlim([55 100])
%ylabel('Probability')
xlabel('Neutralization (%)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Neutralization_2nd_max');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);

%%
Edges_new = 0:5:100;

figure(67)
h1 = histogram(max(hcw_ci_data(:,2001:4000), [], 2),Edges_new, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(max(senior_ci_data(:,2001:4000), [], 2),Edges_new, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
plot([mean(max(hcw_ci_data(:,2001:4000), [], 2)), mean(max(hcw_ci_data(:,2001:4000), [], 2))],[0 h1.Values(20)],'LineWidth',1.5,'Color','#b24a4b');
hold on
plot([mean(max(senior_ci_data(:,2001:4000), [], 2)), mean(max(senior_ci_data(:,2001:4000), [], 2))],[0 h2.Values(20)],'LineWidth',1.5,'Color','#4276bc');
hold on
%ylim([0 0.8])
xlim([55 100])
%ylabel('Probability')
xlabel('Neutralization (%)')
box off
%legend([h2,h1 ],'hcw_ci_data','HHL','Location','Northeast')

figfile = fullfile(pathname,'Neutralization_3rd_max');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI close(gcf);
close(gcf);





