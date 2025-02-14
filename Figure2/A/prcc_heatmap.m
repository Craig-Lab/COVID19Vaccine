clear all
load('Model_LHS.mat');
data = load('prcc_Model_LHS.mat');
prcc = data.prcc;
pvalue = data.sign;
prclab={'peak_A1','peak_A2','peak_time1','peak_time2',' auc_A1','auc_A2','auc_total'};
prctitles={'Antibody peak1','Antibody peak2','Peak time1','Peak time 2','AUC 1','AUC 2',' AUC total'};
PRCC_var={'\delta_{LV}','d_l','d_V','\delta_{TV}','d_T','\delta_{BT}','\rho_G','\rho_S', 'd_B','\delta_{IG}','\beta_{G}',...
    'p_G','d_G','p_{P1}','d_P','d_S','d_M','p_M','\beta_M','p_{P2}','ht','\alpha_P','\alpha_S','d_A', '\rho_I','d_I','SI'};% checked!

pathname=fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
map = {'#678bbe', '#9baecb', '#cfd4e0', '#faf5f4', '#e7c8c6', '#d39794', '#bf6765'};
col_days2 = validatecolor(map, 'multiple');
k= 1;
idx =pvalue{1, k} <0.05;
pp = pvalue{k}(idx);
prcc{k} = prcc{k}(idx);
prcc_var = PRCC_var(idx); 
prcc_val = prcc{k};
[val1, idx1] = sort(abs(prcc_val),'descend'); 
prcc_val1 = prcc_val(idx1);
prcc_var1 = prcc_var(idx1);
[val2, idx2] = sort(prcc_val1,'descend'); 
prcc_val2 = prcc_val1(idx2);
prcc_var2 = prcc_var1(idx2);

[val3, idx3] = sort(abs(prcc_val),'descend'); 
prcc_val3 = prcc_val(idx3(1:10));
prcc_var3 = prcc_var(idx3(1:10));
[val4, idx4] = sort(prcc_val3,'descend'); 
prcc_val4 = prcc_val3(idx4);
prcc_var4 = prcc_var3(idx4);
pp_value_old1 = pp(idx3(1:10));
pp_value1 = pp_value_old1(idx4);
yvalues = {'AP1'};
% Split prcc_var2 and prcc_val2 into two halves
n = ceil(length(prcc_var2) / 2);  % Midpoint index

prcc_var2_top = prcc_var2(1:n);
prcc_val2_top = prcc_val2(1:n);

prcc_var2_bottom = prcc_var2(n+1:end);
prcc_val2_bottom = prcc_val2(n+1:end);

yvalues_top = {'AP1 - Top Half'};
yvalues_bottom = {'AP1 - Bottom Half'};
figure;
t = tiledlayout(2,1, 'TileSpacing', 'none', 'Padding', 'none'); % Ensures uniform spacing
%%
% First Row Heatmap (Top Half)
figure;
t = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
h1 = heatmap(prcc_var2_top, yvalues_top, prcc_val2_top, ...
             'GridVisible', 'off', 'Colormap', col_days2, ...
             'ColorLimits', [-1 1], 'FontColor', 'k');
h1.CellLabelFormat = '%.2f'; % Display values
h1.ColorbarVisible = 'off';
h1.FontSize = 8;
h1.YDisplayLabels = {' '}; % Hide Y-axis labels
h1.XDisplayLabels = prcc_var2_top; % Ensure correct alignment


% Second Row Heatmap (Bottom Half)
nexttile;
h2 = heatmap(prcc_var2_bottom, yvalues_bottom, prcc_val2_bottom, ...
             'GridVisible', 'off', 'Colormap', col_days2, ...
             'ColorLimits', [-1 1], 'FontColor', 'k');
h2.CellLabelFormat = '%.2f'; % Display values
h2.ColorbarVisible = 'off';
h2.FontSize = 8;
h2.YDisplayLabels = {' '}; % Hide Y-axis labels
h2.XDisplayLabels = prcc_var2_bottom; % Ensure correct alignment

% Adjust figure size for two rows
figfile = fullfile(pathname, 'AP1_all_split');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6.5 3]); % Adjusted height for two rows
set(gcf, 'PaperSize', [6.5 3]);

% Save the figure
print(figfile, '-dpng', '-r1080');
close(gcf);
%%
figure;
t = tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
h1 = heatmap(prcc_var2_top, yvalues_top, prcc_val2_top, ...
             'GridVisible', 'off', 'Colormap', col_days2, ...
             'ColorLimits', [-1 1], 'FontColor', 'k');
h1.CellLabelFormat = '%.2f'; % Display values
h1.ColorbarVisible = 'off';
h1.FontSize = 8;
h1.YDisplayLabels = {' '}; % Hide Y-axis labels
h1.XDisplayLabels = prcc_var2_top; % Ensure correct alignment

% Adjust figure size for two rows
figfile = fullfile(pathname, 'AP1_1');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6.5 1.5]); % Adjusted height for two rows
set(gcf, 'PaperSize', [6.5 1.5]);
print(figfile, '-dpng', '-r1080');
close(gcf);

figure;
t = tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
h2 = heatmap(prcc_var2_bottom, yvalues_bottom, prcc_val2_bottom, ...
             'GridVisible', 'off', 'Colormap', col_days2, ...
             'ColorLimits', [-1 1], 'FontColor', 'k');
h2.CellLabelFormat = '%.2f'; % Display values
h2.ColorbarVisible = 'off';
h2.FontSize = 8;
h2.YDisplayLabels = {' '}; % Hide Y-axis labels
h2.XDisplayLabels = prcc_var2_bottom; % Ensure correct alignment

% Adjust figure size for two rows
figfile = fullfile(pathname, 'AP1_2');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6.5 1.5]); % Adjusted height for two rows
set(gcf, 'PaperSize', [6.5 1.5]);

% Save the figure
print(figfile, '-dpng', '-r1080');
close(gcf);

%%
h = heatmap(prcc_var2,yvalues,prcc_val2, 'GridVisible','off', 'Colormap',col_days2,'ColorLimits',[-1 1],'FontColor','k');
h.CellLabelFormat = '%.2f';
h.ColorbarVisible = 'off';
h.CellLabelColor = 'k';
h.FontSize =8;
h.InnerPosition(4) = h.InnerPosition(4) - 0.4; % Reduce height slightly
h.Position(2) = h.Position(2) + 0.4;          % Push heatmap upward
h.YDisplayLabels = {' '};
h.GridVisible = 'off';
figfile = fullfile(pathname,'AP1_all');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0  22 2]); 
set(gcf, 'PaperSize', [22 2]);
print(figfile, '-dpng', '-r1080');
close(gcf);
%%
h = heatmap(prcc_var4,yvalues,prcc_val4, 'GridVisible','off', 'Colormap',col_days2,'ColorLimits',[-1 1],'FontColor','k');
h.CellLabelFormat = '%.2f';
h.ColorbarVisible = 'off';
h.CellLabelColor = 'k';
h.FontSize =8;
h.InnerPosition(4) = h.InnerPosition(4) - 0.4; % Reduce height slightly
h.Position(2) = h.Position(2) + 0.4;          % Push heatmap upward
h.YDisplayLabels = {' '};
h.GridVisible = 'off';
figfile = fullfile(pathname,'AP1');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0  12 2]); 
set(gcf, 'PaperSize', [12 2]);
print(figfile, '-dpng', '-r1080');
close(gcf);

k = 2;
idx = pvalue{1, k} < 0.01;
pp = pvalue{k}(idx);
prcc{k} = prcc{k}(idx);
prcc_var = PRCC_var(idx);
prcc_val = prcc{k};

% Sort by absolute PRCC values in descending order
[val1, idx1] = sort(abs(prcc_val), 'descend');
prcc_val1 = prcc_val(idx1);
prcc_var1 = prcc_var(idx1);

% Sort the top values in descending order
[val2, idx2] = sort(prcc_val1, 'descend');
prcc_val2 = prcc_val1(idx2);
prcc_var2 = prcc_var1(idx2);

% Select top 10 absolute PRCC values
[val3, idx3] = sort(abs(prcc_val), 'descend');
prcc_val3 = prcc_val(idx3(1:10));
prcc_var3 = prcc_var(idx3(1:10));

[val4, idx4] = sort(prcc_val3, 'descend');
prcc_val4 = prcc_val3(idx4);
prcc_var4 = prcc_var3(idx4);
pp_value_old2 = pp(idx3(1:10));
pp_value2 = pp_value_old2(idx4);

% Split the data into two halves for two rows
split_idx = ceil(length(prcc_var2) / 2);
prcc_var_top = prcc_var2(1:split_idx);
prcc_val_top = prcc_val2(1:split_idx);
prcc_var_bottom = prcc_var2(split_idx+1:end);
prcc_val_bottom = prcc_val2(split_idx+1:end);
%
%%

figure;
t = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top Heatmap
nexttile;
h1 = heatmap(prcc_var_top, yvalues, prcc_val_top, 'GridVisible', 'off', ...
    'Colormap', col_days2, 'ColorLimits', [-1 1], 'FontColor', 'k');
h1.CellLabelFormat = '%.2f';
h1.ColorbarVisible = 'off';
h1.FontSize = 8;
h1.YDisplayLabels = {' '}; % Set to appropriate labels if needed
h1.XDisplayLabels = prcc_var_top; 

nexttile;
h2 = heatmap(prcc_var_bottom, yvalues, prcc_val_bottom, 'GridVisible', 'off', ...
    'Colormap', col_days2, 'ColorLimits', [-1 1], 'FontColor', 'k');
h2.CellLabelFormat = '%.2f';
h2.ColorbarVisible = 'off';
h2.FontSize = 8;
h2.YDisplayLabels = {' '}; % Set to appropriate labels if needed
h2.XDisplayLabels = prcc_var_bottom;


% Save the figure
figfile = fullfile(pathname, 'AP2_all');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 9 3]); 
set(gcf, 'PaperSize', [9 3]);
print(figfile, '-dpng', '-r1080');

close(gcf);
%%


figure;
t = tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
h1 = heatmap(prcc_var_top, yvalues, prcc_val_top, 'GridVisible', 'off', ...
    'Colormap', col_days2, 'ColorLimits', [-1 1], 'FontColor', 'k');
h1.CellLabelFormat = '%.2f';
h1.ColorbarVisible = 'off';
h1.FontSize = 8;
h1.YDisplayLabels = {' '}; % Set to appropriate labels if needed
h1.XDisplayLabels = prcc_var_top; 

figfile = fullfile(pathname, 'AP2_1');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 9 1.5]); 
set(gcf, 'PaperSize', [9 1.5]);
print(figfile, '-dpng', '-r1080');
close(gcf);
%%
figure;
t = tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
h2 = heatmap(prcc_var_bottom, yvalues, prcc_val_bottom, 'GridVisible', 'off', ...
    'Colormap', col_days2, 'ColorLimits', [-1 1], 'FontColor', 'k');

h2.CellLabelFormat = '%.2f'; % Force numbers to show
h2.CellLabelColor = 'k'; % Ensure text is visible
h2.ColorbarVisible = 'off';
h2.FontSize = 8; % Increase font size
h2.YDisplayLabels = {' '}; 
h2.XDisplayLabels = prcc_var_bottom; 

% Save the figure
figfile = fullfile(pathname, 'AP2_2');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8.3 1.5]); 
set(gcf, 'PaperSize', [8.3 1.5]);
print(figfile, '-dpng', '-r1080');

close(gcf);




%%
k= 2;
idx =pvalue{1, k} <0.01;
pp = pvalue{k}(idx);
prcc{k} = prcc{k}(idx);
prcc_var = PRCC_var(idx); 
prcc_val = prcc{k};
[val1, idx1] = sort(abs(prcc_val),'descend'); 
prcc_val1 = prcc_val(idx1);
prcc_var1 = prcc_var(idx1);
[val2, idx2] = sort(prcc_val1,'descend'); 
prcc_val2 = prcc_val1(idx2);
prcc_var2 = prcc_var1(idx2);

[val3, idx3] = sort(abs(prcc_val),'descend'); 
prcc_val3 = prcc_val(idx3(1:10));
prcc_var3 = prcc_var(idx3(1:10));

[val4, idx4] = sort(prcc_val3,'descend'); 
prcc_val4 = prcc_val3(idx4);
prcc_var4 = prcc_var3(idx4);
pp_value_old2 = pp(idx3(1:10));
pp_value2 = pp_value_old2(idx4);
yvalues = {'AP2'};
h = heatmap(prcc_var2,yvalues,prcc_val2, 'GridVisible','off', 'Colormap',col_days2,'ColorLimits',[-1 1],'FontColor','k');
h.CellLabelFormat = '%.2f';
h.ColorbarVisible = 'off';
h.CellLabelColor = 'k';
h.FontSize =8;
h.InnerPosition(4) = h.InnerPosition(4) - 0.4; % Reduce height slightly
h.Position(2) = h.Position(2) + 0.4;          % Push heatmap upward
h.YDisplayLabels = {' '};
h.GridVisible = 'off';
figfile = fullfile(pathname,'AP2_all');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0  20 2]); 
set(gcf, 'PaperSize', [20 2]);
print(figfile, '-dpng', '-r1080');
close(gcf);
%%
h = heatmap(prcc_var4,yvalues,prcc_val4, 'GridVisible','off', 'Colormap',col_days2,'ColorLimits',[-1 1],'FontColor','k');
h.CellLabelFormat = '%.2f';
h.ColorbarVisible = 'off';
h.CellLabelColor = 'k';
h.FontSize =8;
h.InnerPosition(4) = h.InnerPosition(4) - 0.4; % Reduce height slightly
h.Position(2) = h.Position(2) + 0.4;          % Push heatmap upward
h.YDisplayLabels = {' '};
h.GridVisible = 'off';
figfile = fullfile(pathname,'AP2');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0  12 2]); 
set(gcf, 'PaperSize', [12 2]);
print(figfile, '-dpng', '-r1080');
close(gcf);
%%
close all
load('Model_LHS.mat');
data = load('prcc_Model_LHS.mat');
prcc = data.prcc;
pvalue = data.sign;
prclab={'peak_A1','peak_A2','peak_time1','peak_time2',' auc_A1','auc_A2','auc_total'};
prctitles={'Antibody peak1','Antibody peak2','Peak time1','Peak time 2','AUC 1','AUC 2',' AUC total'};
PRCC_var={'\delta_{LV}','d_l','d_V','\delta_{TV}','d_T','\delta_{BT}','\rho_G','\rho_S', 'd_B','\delta_{IG}','\beta_{G}',...
    'p_G','d_G','p_{P1}','d_P','d_S','d_M','p_M','\beta_M','p_{P2}','ht','\alpha_P','\alpha_S','d_A', '\rho_I','d_I','SI'};% checked!
corr_value = 0.7;
PRCC_PLOT(LHSmatrix, peak_A1, 1, PRCC_var, 'Antibody peak', corr_value, 'peak_A1')
close all
%%
close all
load('Model_LHS.mat');
data = load('prcc_Model_LHS.mat');
prcc = data.prcc;
pvalue = data.sign;
prclab={'peak_A1','peak_A2','peak_time1','peak_time2',' auc_A1','auc_A2','auc_total'};
prctitles={'Antibody peak1','Antibody peak2','Peak time1','Peak time 2','AUC 1','AUC 2',' AUC total'};
PRCC_var={'\delta_{LV}','d_l','d_V','\delta_{TV}','d_T','\delta_{BT}','\rho_G','\rho_S', 'd_B','\delta_{IG}','\beta_{G}',...
    'p_G','d_G','p_{P1}','d_P','d_S','d_M','p_M','\beta_M','p_{P2}','ht','\alpha_P','\alpha_S','d_A', '\rho_I','d_I','SI'};% checked!
corr_value = 0.29;
PRCC_PLOT(LHSmatrix, peak_A2, 1, PRCC_var, 'Antibody peak', corr_value,'peak_A2')
close all
