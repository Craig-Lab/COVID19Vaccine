%%
clear all
format long
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
threshold = [60, 70, 80, 90, 95];
%% timing for the fourth dose (half third dose)
load('neu_senior_LHH.mat')
load('neu_senior_HHH.mat')




I1= zeros(27,5); % HHH_neu
I2= zeros(27,5); % LHH_neu


[a1, IP1] = max(HHH_neu(:,2001:4000),[],2); %HHH_neu
[a2, IP2] = max(LHH_neu(:,2001:4000),[],2); %LHH_neu

IP1 = IP1+2000;
IP2 = IP2+2000;

%
for i = 1:27
    for j = 1:5
            idx = IP1(i);
                if HHH_neu(i,idx)<threshold(j)
                    I1(i,j) =2001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(HHH_neu(i,idx:4000)>threshold(j)));
                     I1(i,j)= a+idx-1;
                end
    end
end
for i = 1:27
    for j = 1:5
            idx = IP2(i);
                if LHH_neu(i,idx)<threshold(j)
                    I2(i,j) =2001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(LHH_neu(i,idx:4000)>threshold(j)));
                     I2(i,j)= a+idx-1;
                end
    end
end


T1 = round(time(I1));
T2 = round(time(I2));

T1(T1 >= 629)= 630;
T2(T2 >= 629)= 630;

T1 = T1-270;
T2 = T2-270;

% compare median of those who require additional dose within one year
m1 =median(T1);
m2 = median(T2);

m = [m1,median(HHH_neu(:,4000));m2,median(LHH_neu(:,4000))];

%%
figure(31)
h1 = histogram(HHH_neu(:,4000),0:5:100, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on 
h2 = histogram(LHH_neu(:,4000),0:5:100, 'Normalization','Probability','FaceColor','#3E7748','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(HHH_neu(:,4000)), median(HHH_neu(:,4000))],[0 h1.Values(14)+0.02],':','LineWidth',1.5, 'Color','#2e4a85');
hold on
plot([median(LHH_neu(:,4000)), median(LHH_neu(:,4000))],[0 h2.Values(14)],'LineWidth',1.5,'Color',[0.24, 0.47, 0.28, 1]);
hold on
ylim([0 0.3])
xlim([30 100])
yticks([0, 0.1, 0.2, 0.3])
%ylabel('Probability')
xlabel('Neutralization (%)')
box off
figfile = fullfile(pathname,'LHH_neu_dist_senior');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 

set(gcf, 'PaperPosition', [0 0 3.8 3.2]); 
set(gcf, 'PaperSize', [3.8 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);
%%
Edges= 0:20:380;
figure(31)
h1 = histogram(T1(:,1),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on 
h2 = histogram(T2(:,1),Edges, 'Normalization','Probability','FaceColor','#3E7748','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m2(1)), median(m2(1))],[0 h2.Values(19)],'LineWidth',1.5, 'Color',[0.24, 0.47, 0.28, 0.9]);
hold on
plot([median(m1(1)), median(m1(1))],[0 h1.Values(19)+0.02],':','LineWidth',1.5, 'Color','#2e4a85');
hold on
ylim([0 0.6])
xlim([0 380])
%yticks([0,0.1,0.2,0.5])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LHH_neu_60_senior');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3 3]); 
set(gcf, 'PaperSize', [ 3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);
%%
figure(32)
h1 = histogram(T1(:,2),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on 
h2 = histogram(T2(:,2),Edges, 'Normalization','Probability','FaceColor','#3E7748','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m2(2)), median(m2(2))],[0 h1.Values(17)],'LineWidth',1.5, 'Color',[0.24, 0.47, 0.28,1]);
hold on
plot([median(m1(2)), median(m1(2))],[0 0.2],':','LineWidth',1, 'Color','#2e4a85');
hold on

ylim([0 0.5])
xlim([0 380])
yticks([0,0.1,0.2,0.5])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LHH_neu_70_senior');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0  3 3]); 
set(gcf, 'PaperSize', [ 3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);

%%
%
figure(33)
h1 = histogram(T1(:,3),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(T2(:,3),Edges, 'Normalization','Probability','FaceColor','#3E7748','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m2(3)), median(m2(3))],[0 h1.Values(14)],'LineWidth',1.5, 'Color',[0.24, 0.47, 0.28,1]);
hold on
plot([median(m1(3)), median(m1(3))],[0 0.15],':','LineWidth',1, 'Color','#2e4a85');

ylim([0 0.4])
yticks(0:0.2:0.4)
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_LHH_neu_80_senior');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0  3 3]); 
set(gcf, 'PaperSize', [ 3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%
figure(34)
h1 = histogram(T1(:,4),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.3);
hold on
h2 = histogram(T2(:,4),Edges, 'Normalization','Probability','FaceColor','#3E7748','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m2(4)), median(m2(4))],[0 h1.Values(15)],'LineWidth',1.5, 'Color',[0.24, 0.47, 0.28,1]);
hold on
plot([median(m1(4)), median(m1(4))],[0 0.15],':','LineWidth',1, 'Color','#2e4a85');
hold on
xlim([0 380])
ylim([0 0.3])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LHH_neu_90_senior');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0  3 3]); 
set(gcf, 'PaperSize', [ 3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

%%
figure(35)
h1 = histogram(T1(:,5),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.25);
hold on
h2 = histogram(T2(:,5),Edges, 'Normalization','Probability','FaceColor','#3E7748','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m2(5)), median(m2(5))],[0 h1.Values(7)],'LineWidth',1.5, 'Color',[0.24, 0.47, 0.28,1]);
hold on
plot([median(m1(5)), median(m1(5))],[0 0.15],':','LineWidth',1, 'Color','#2e4a85');
hold on
xlim([0 380])
ylim([0 0.3])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LHH_neu_95_senior');
set(gca,'Fontsize',7)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0  3 3]); 
set(gcf, 'PaperSize', [ 3 3]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);

