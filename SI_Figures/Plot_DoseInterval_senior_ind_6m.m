clear all
format long
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
SH= load('Neu_senior_int12_6m.mat').SH;
SA_data= load('Dose_6m_senior.mat').SA;
time  =  load('Neu_senior_int12_6m.mat').time;


%%
threshold = [60, 70, 80, 90, 95];

% timing for the fourth dose 
SSdata=load('senior_neu_ind.mat');
SS = SSdata.SS;
I1= zeros(27,5); % 2m
I2= zeros(27,5); % 6m

IP1= zeros(27,1); %2m
IP2= zeros(27,1); %6m
for i= 1:27
    [a1, b1] = max(SS(i,2001:4000)); %2m
    [a2, b2] = max(SH(i,2001:4000)); %6m
    IP1(i) = b1+2000;
    IP2(i) = b2+2000;    
end

for i = 1:27
    for j = 1:5
            idx = IP1(i);
                if SS(i,idx)<threshold(j)
                    I1(i,j) =2001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(SS(i,idx:4000)>threshold(j)));
                     I1(i,j)= a+idx-1;
                end
    end
end
for i = 1:27
    for j = 1:5
            idx = IP2(i);
                if SH(i,idx)<threshold(j)
                    I2(i,j) =2001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(SH(i,idx:4000)>threshold(j)));
                     I2(i,j)= a+idx-1;
                end
    end
end


T1 = round(time(I1));
T2 = round(time(I2));

T1(T1 >= 629)= 630;
T2(T2 >= 629)= 630;

%T1 = T1-270;
%T2 = T2-270;

% compare median of those who require additional dose within one year
m1 =median(T1);
m2 = median(T2);

m = [m1;m2]; 

a1 = sort(T1(:,4));
a2 = sort(T2(:,4));


figure(31)
h1 = histogram(SS(:,4000),0:10:100, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.2);
hold on 
h2 = histogram(SH(:,4000),0:10:100, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(SS(:,4000)), median(SS(:,4000))],[0 h1.Values(7)],'LineWidth',1.5, 'Color','#4276bc');
hold on
plot([median(SH(:,4000)), median(SH(:,4000))],[0 h2.Values(9)],'LineWidth',1.5, 'Color','#3e7749');
ylim([0 0.4])
xlim([30 100])
yticks([0, 0.2, 0.4])
%xlabel('Time (days)')
xlabel('Neutralization (%)')
box off
figfile = fullfile(pathname,'6m_dist_senior');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.5 3.2]); 
set(gcf, 'PaperSize', [3.5 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);

%%
Edges= 0:20:700;

figure(31)
h1 = histogram(T1(:,1),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.2);
hold on 
h2 = histogram(T2(:,1),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m1(1)), median(m1(1))],[0 h1.Values(32)],'LineWidth',1.5, 'Color','#4276bc');
hold on
plot([median(m2(1)), median(m2(1))],[0 h2.Values(32)],'LineWidth',1.5, 'Color','#3e7749');
ylim([0 1])
xlim([400 650])
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_6m_60');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.5 3.2]); 
set(gcf, 'PaperSize', [3.5 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);
%%
figure(32)
h1 = histogram(T1(:,2),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.2);
hold on 
h2 = histogram(T2(:,2),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m1(2)), median(m1(2))],[0 h1.Values(31)],'LineWidth',1.5, 'Color','#4276bc');
hold on
plot([median(m2(2)), median(m2(2))],[0 h2.Values(32)],'LineWidth',1.5, 'Color','#3e7749');
ylim([0 0.8])
yticks(0:0.4:0.8)
xlim([400 650])
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_6m_70');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.5 3.2]); 
set(gcf, 'PaperSize', [3.5 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);

%%
figure(33)
h1 = histogram(T1(:,3),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.2)
hold on
h2 = histogram(T2(:,3),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.25)
hold on
plot([median(m1(3)), median(m1(3))],[0 h1.Values(29)],'LineWidth',1.5, 'Color','#4276bc');
hold on
plot([median(m2(3)), median(m2(3))],[0 h2.Values(32)],'LineWidth',1.5, 'Color','#3e7749');
hold off
ylim([0 0.6])
xlim([400 650])
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_6m_80');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.5 3.2]); 
set(gcf, 'PaperSize', [3.5 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);

%%
figure(34)
h1 = histogram(T1(:,4),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.2);
hold on
h2 = histogram(T2(:,4),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m1(4)), median(m1(4))],[0 h1.Values(26)],'LineWidth',1.5, 'Color','#4276bc');
hold on
plot([median(m2(4)), median(m2(4))],[0 h2.Values(28)],'LineWidth',1.5, 'Color','#3e7749');
hold on
ylim([0 0.4])
xlim([400 650])
%ylabel('Probability')
xlabel('Time (days)')
box off
set(gca,'Fontsize',12)
figfile = fullfile(pathname,'Time3to4_6m_90');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.5 3.2]); 
set(gcf, 'PaperSize', [3.5 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);

%
%%
figure(35)
h1 = histogram(T1(:,5),Edges, 'Normalization','Probability','FaceColor','#4276bc','EdgeAlpha',0,'facealpha', 0.2);
hold on
h2 = histogram(T2(:,5),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.25);
hold on
plot([median(m1(5)), median(m1(5))],[0 h1.Values(24)],'LineWidth',1.5, 'Color','#4276bc');
hold on
plot([median(m2(5)), median(m2(5))],[0 h2.Values(26)],'LineWidth',1.5, 'Color','#3e7749');
hold on
ylim([0 0.2])
yticks(0:0.1:0.2)
xlim([400 650])
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_6m_95');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gcf, 'PaperPosition', [0 0 3.5 3.2]); 
set(gcf, 'PaperSize', [3.5 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r800');  % Save as PNG at 300 DPI
close(gcf);
%%
a = sort(T1(:,4),'descend');
b = sort(T2(:,4),'descend');
c = sort(SH(:,4000),'descend');
sum(a==630)
sum(b==630)
%%
median(SS(:,4000))
median(SH(:,4000))