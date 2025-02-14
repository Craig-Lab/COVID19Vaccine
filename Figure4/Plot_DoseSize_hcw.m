%%
clear all
format long
pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');

threshold = [60, 70, 80, 90, 95];
% timing for the fourth dose (half third dose)
load('neu_HCW_LHH.mat')
load('neu_HCW_LLH.mat')
load('neu_HCW_HHL.mat')



I1= zeros(32,5); % HHL_neu
I2= zeros(32,5); % LHH_neu
I3= zeros(32,5); % LLH_neu

HHL_neu_ci = zeros(2,4000);
LHH_neu_ci = zeros(2,4000);
LLH_neu_ci = zeros(2,4000); 
a_end = [median(HHL_neu(:,4000)),median(LHH_neu(:,4000)),median(LLH_neu(:,4000))];
a_6m = [median(HHL_neu(:,3001)),median(LHH_neu(:,3001)),median(LLH_neu(:,3001))];

[a1, IP1] = max(HHL_neu(:,2001:4000),[],2); %HHL_neu
[a2, IP2] = max(LHH_neu(:,2001:4000),[],2); %LHH_neu
[a3, IP3] = max(LLH_neu(:,2001:4000),[],2);% LLH_neu
IP1 = IP1+2000;
IP2 = IP2+2000;
IP3 = IP3+2000;
%
for i = 1:32
    for j = 1:5
            idx = IP1(i);
                if HHL_neu(i,idx)<threshold(j)
                    I1(i,j) =2001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(HHL_neu(i,idx:4000)>threshold(j)));
                     I1(i,j)= a+idx-1;
                end
    end
end
for i = 1:32
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

for i = 1:32
    for j = 1:5
            idx = IP3(i);
                if LLH_neu(i,idx)<threshold(j)
                    I3(i,j) =2001; %% peak value cannot reach the corresponding threshold
                else
                    a = max(find(LLH_neu(i,idx:4000)>threshold(j)));
                     I3(i,j)= a+idx-1;
                end
    end
end
T1 = round(time(I1));
T2 = round(time(I2));
T3 = round(time(I3));
T1(T1 >= 629)= 630;
T2(T2 >= 629)= 630;
T3(T3 >= 629)= 630;
T1 = T1-270;
T2 = T2-270;
T3 = T3-270;
% compare median of those who require additional dose within one year
m1 =median(T1);
m2 = median(T2);
m3 = median(T3);
m = [m1;m2;m3];

%%
a = T1(:,2:5);
b = T2(:,2:5);%%LHH
c = T3(:,2:5);%%LLH
[~, idx] = sort(a(:, 1), 'descend');
a_new = a(idx,:);
b_new = b(idx,:);
c_new = c(idx,:);
aa = sort(a,'descend');
bb = sort(b,'descend');
cc = sort(c,'descend');
d1 =b_new-a_new;
d2= c_new-a_new;
%%

figure(111)
h1 = histogram(HHL_neu(:,4000),0:5:100, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on 
h2 = histogram(LHH_neu(:,4000),0:5:100, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(a_end(1)), median(a_end(1))],[0 h1.Values(13)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(a_end(2)), median(a_end(2))],[0 h2.Values(17)],'LineWidth',1.5, 'Color','#3e7749');
ylim([0 0.3])
xlim([30 100])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'LHH_neu_dist_end');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
% set(gcf, 'PaperPosition', [0 0  3 3]); 
% set(gcf, 'PaperSize', [ 3 3]);

set(gcf, 'PaperPosition', [0 0 3 3.2]); 
set(gcf, 'PaperSize', [3 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);
%%
figure(113)
h1 = histogram(HHL_neu(:,4000),0:5:100, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on 
h2 = histogram(LLH_neu(:,4000),0:5:100, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(a_end(1)), median(a_end(1))],[0 h1.Values(13)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(a_end(3)), median(a_end(3))],[0 h2.Values(16)],'LineWidth',1.5, 'Color','#3e7749');
ylim([0 0.3])
xlim([30 100])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Neutralization(%)')
box off
figfile = fullfile(pathname,'LLH_neu_dist_end');
set(gca,'Fontsize',8)
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
%set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
% set(gcf, 'PaperPosition', [0 0  3 3]); 
% set(gcf, 'PaperSize', [ 3 3]);

set(gcf, 'PaperPosition', [0 0 3 3.2]); 
set(gcf, 'PaperSize', [3 3.2]);
set(gca, 'LooseInset', get(gca,'TightInset'))
saveas(gcf,figfile,'pdf')
print(figfile, '-dpng', '-r1080');  % Save as PNG at 300 DPI
close(gcf);


%%
Edges= 0:20:380;
figure(32)
h1 = histogram(T1(:,1),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on 
h2 = histogram(T2(:,1),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(1)), median(m1(1))],[0 h1.Values(19)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m2(1)), median(m2(1))],[0 h2.Values(19)],'LineWidth',1.5, 'Color','#3e7749');
%ylim([0 0.6])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LHH_neu_60');
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
figure(32)
h1 = histogram(T1(:,2),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on 
h2 = histogram(T2(:,2),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(2)), median(m1(2))],[0 h1.Values(17)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m2(2)), median(m2(2))],[0 h2.Values(19)],'LineWidth',1.5, 'Color','#3e7749');
%ylim([0 0.6])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LHH_neu_70');
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
figure(33)
h1 = histogram(T1(:,3),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.25)
hold on
h2 = histogram(T2(:,3),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(3)), median(m1(3))],[0 h1.Values(13)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m2(3)), median(m2(3))],[0 h2.Values(19)],'LineWidth',1.5, 'Color','#3e7749');
hold off
ylim([0 0.6])
%yticks([0,0.1,0.2,0.5])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_LHH_neu_80');
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
h1 = histogram(T1(:,4),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on
h2 = histogram(T2(:,4),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(4)), median(m1(4))],[0 h1.Values(10)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m2(4)), median(m2(4))],[0 h2.Values(15)],'LineWidth',1.5, 'Color','#3e7749');
hold on
xlim([0 380])
ylim([0 0.3])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Time (days)')
box off
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LHH_neu_90');
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
h1 = histogram(T1(:,5),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on
h2 = histogram(T2(:,5),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(5)), median(m1(5))],[0 h1.Values(7)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m2(5)), median(m2(5))],[0 h2.Values(11)],'LineWidth',1.5, 'Color','#3e7749');
hold on
xlim([0 380])
ylim([0 0.3])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Time (days)')
box off
%legend([h2,h1 ],'LHH_neu','HHL_neu','Location','Northeast')

figfile = fullfile(pathname,'Time3to4_LHH_neu_95');
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
Edges= 0:20:380;
figure(42)
h1 = histogram(T1(:,1),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(T3(:,1),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(2)), median(m1(1))],[0 h1.Values(17)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m3(1)), median(m3(1))],[0 h2.Values(19)],'LineWidth',1.5, 'Color','#3e7749');
%ylim([0 0.6])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LLH_neu_60');
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
figure(42)
h1 = histogram(T1(:,2),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(T3(:,2),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(2)), median(m1(2))],[0 h1.Values(17)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m3(2)), median(m3(2))],[0 h2.Values(19)],'LineWidth',1.5, 'Color','#3e7749');
ylim([0 0.6])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LLH_neu_70');
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

figure(43)
h1 = histogram(T1(:,3),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.25)
hold on
h2 = histogram(T3(:,3),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(3)), median(m1(3))],[0 h1.Values(13)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m3(3)), median(m3(3))],[0 h2.Values(18)],'LineWidth',1.5, 'Color','#3e7749');
hold off
ylim([0 0.5])
yticks([0 0.1 0.2 0.5])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_LLH_neu_80');
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
figure(44)
h1 = histogram(T1(:,4),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on
h2 = histogram(T3(:,4),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(4)), median(m1(4))],[0 h1.Values(10)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m3(4)), median(m3(4))],[0 h2.Values(14)],'LineWidth',1.5, 'Color','#3e7749');
hold on
xlim([0 380])
ylim([0 0.3])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Time (days)')
box off

figfile = fullfile(pathname,'Time3to4_LLH_neu_90');
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
figure(45)
h1 = histogram(T1(:,5),Edges, 'Normalization','Probability','FaceColor','#b24a4b','EdgeAlpha',0,'facealpha', 0.28);
hold on
h2 = histogram(T3(:,5),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m1(5)), median(m1(5))],[0 h1.Values(7)],'LineWidth',1.5, 'Color','#b24a4b');
hold on
plot([median(m3(5)), median(m3(5))],[0 h2.Values(10)],'LineWidth',1.5, 'Color','#3e7749');
hold on
xlabel('Time (days)')
box off
xlim([0 380])
ylim([0 0.3])
yticks(0:0.1:0.3)
figfile = fullfile(pathname,'Time3to4_LLH_neu_95');
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
Edges= 0:20:380;
figure(52)
h1 = histogram(T2(:,1),Edges, 'Normalization','Probability','FaceColor','#375391','EdgeAlpha',0,'facealpha', 0.28);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(T3(:,1),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m2(1)), median(m2(1))],[0 h1.Values(19)],'LineWidth',1.5, 'Color','#375391');
hold on
plot([median(m3(1)), median(m3(1))],[0 h2.Values(19)],'LineWidth',1.5, 'Color','#3e7749');
%ylim([0 0.6])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LLH_LHH_neu_60');
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
figure(52)
h1 = histogram(T2(:,2),Edges, 'Normalization','Probability','FaceColor','#375391','EdgeAlpha',0,'facealpha', 0.28);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(T3(:,2),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m2(2)), median(m2(2))],[0 h1.Values(19)],'LineWidth',1.5, 'Color','#375391');
hold on
plot([median(m3(2)), median(m3(2))],[0 h2.Values(19)],'LineWidth',1.5, 'Color','#3e7749');
%ylim([0 0.6])
xlim([0 380])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LLH_LHH_neu_70');
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
figure(53)
h1 = histogram(T2(:,3),Edges, 'Normalization','Probability','FaceColor','#375391','EdgeAlpha',0,'facealpha', 0.28);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(T3(:,3),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m2(3)), median(m2(3))],[0 h1.Values(19)],'LineWidth',1.5, 'Color','#375391');
hold on
plot([median(m3(3)), median(m3(3))],[0 h2.Values(18)],'LineWidth',1.5, 'Color','#3e7749');
ylim([0 0.5])
xlim([0 380])
yticks([0 0.1 0.2 0.5])
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LLH_LHH_neu_80');
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

figure(54)
h1 = histogram(T2(:,4),Edges, 'Normalization','Probability','FaceColor','#375391','EdgeAlpha',0,'facealpha', 0.28);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(T3(:,4),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m2(4)), median(m2(4))],[0 h1.Values(15)],'LineWidth',1.5, 'Color','#375391');
hold on
plot([median(m3(4)), median(m3(4))],[0 h2.Values(14)],'LineWidth',1.5, 'Color','#3e7749');
hold off
xlim([0 380])
ylim([0 0.3])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LLH_LHH_neu_90');
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

figure(55)
h1 = histogram(T2(:,5),Edges, 'Normalization','Probability','FaceColor','#375391','EdgeAlpha',0,'facealpha', 0.28);
hold on %#b24a4b');%'#b24a4b'
%text(10, 30, ['Median = ' num3str(round(median(a1),2))], 'FontSize', 14)
h2 = histogram(T3(:,5),Edges, 'Normalization','Probability','FaceColor','#3e7749','EdgeAlpha',0,'facealpha', 0.28);
hold on
plot([median(m2(5)), median(m2(5))],[0 h1.Values(11)],'LineWidth',1.5, 'Color','#375391');
hold on
plot([median(m3(5)), median(m3(5))],[0 h2.Values(10)],'LineWidth',1.5, 'Color','#3e7749');
hold off
xlim([0 380])
ylim([0 0.3])
yticks(0:0.1:0.3)
%ylabel('Probability')
xlabel('Time (days)')
box off
figfile = fullfile(pathname,'Time3to4_LLH_LHH_neu_95');
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
median(LHH_neu(:,4000))
median(LLH_neu(:,4000))
median(HHL_neu(:,4000))