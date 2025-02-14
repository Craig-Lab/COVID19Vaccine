%% Plot the residual of the partial regression of X (input: LHS matrix) and Y (output)
%% at column s (one time point saved). RCC Coefficients are calculated on these
%% var: labels of the parameters varied in the X (as legend)
%% The Title of the plot is the Pearson correlation coefficient of the
%% transformed data, that is  the PRCC calculated on the original data.
%% The p-value is also showed in the title
%% By Simeone Marino, June 5 2007 %%
%% Modified by Marissa Renardy and Paul Wolberg, June 8, 2022.
%
% Inputs:
%    X: The LHS Matrix, N x k, where N is the number of runs and k is the
%       number of varied parameters.
%
%   Y: The model outputs, T x N, where T is the number of time points and
%      N is the number of runs.
%
%  s: A single time point ordinal. If T is the number of time points then s is
%     a single value in the range [1, T], i.e. 1 <= s <= T. For example if T is
%     10 then s is in the range [1, 10], i.e. 1 <= s <= 10.
%
% PRCC_var: A cell array of string names of the k varied parameters. This is
%           from the settings file, and is in the Matlab workspace and result
%           .mat file (Model_LHS.mat) after running Model_LHS.m.
%
% y_var: A cell array of string names of the model outputs. This is from the
%        settings file,  and is in the Matlab workspace and result
%        .mat file (Model_LHS.mat) after running Model_LHS.m, with name
%        y_var_label.
%
% For example:
% PRCC_PLOT(LHSmatrix, V_lhs, 1, PRCC_var, y_var_label)
% PRCC_PLOT(LHSmatrix, V_lhs, 2, PRCC_var, y_var_label)

function PRCC_PLOT(X, Y, s, PRCC_var, y_var, corr_value, metric)

Y=Y(s,:);
[a k]=size(X); % Define the size of LHS matrix
Xranked=rankingN(X);
Yranked=ranking1(Y);
for i=1:k  % Loop for the whole submatrices, Zi
    c1=['LHStemp=Xranked;LHStemp(:,',num2str(i),')=[];Z',num2str(i),'=[ones(a,1) LHStemp];LHStemp=[];'];
    eval(c1);
end
for i=1:k
    c2=['[b',num2str(i),',bint',num2str(i),',r',num2str(i),']= regress(Yranked,Z',num2str(i),');'];
    c3=['[b',num2str(i),',bint',num2str(i),',rx',num2str(i),']= regress(Xranked(:,',num2str(i),'),Z',num2str(i),');'];
    eval(c2);
    eval(c3);
end
for i=1:k
    c4=['r',num2str(i)];
    c5=['rx',num2str(i)];
    [r p]=corr(eval(c4),eval(c5));
    %a=['[PRCC , p-value] = ' '[' num2str(r) ' , '  num2str(round(p,2)) '].'];% ' Time point=' num2str(s-1)];

    %a=['PRCC =  ' num2str(round(r,2)) ' ,  P-value \approx '  num2str(round(p,2)) ];% ' Time point=' num2str(s-1)];    
    mdl = fitlm((eval(c4)),(eval(c5)));

    % Get slope and intercept from the linear model
    coeffs = mdl.Coefficients.Estimate;
    slope = coeffs(2);  % The slope (coefficient of x)
    intercept = coeffs(1);  % The intercept
    
    % Define x-limits for the line to restrict it within the data range
    x_fit = linspace(min(eval(c4)), max(eval(c4)), 1000);
    y_fit = slope * x_fit + intercept;
    % Restrict the y-values and corresponding x-values to be within -500 and 500
    valid_idx = (y_fit >= -700) & (y_fit <= 650);
    x_fit = x_fit(valid_idx);
    y_fit = y_fit(valid_idx);
    % Check conditions for plotting based on correlation and p-value
    if p < 0.05 && (r < -corr_value || r > corr_value)
        % Create figure
        fig(i) = figure;
        
        % Select color based on positive or negative correlation
        if r < -corr_value
            color = '#7896C1';  % Color for negative correlation
        else
            color = '#D39794';  % Color for positive correlation
        end
        
        % Plot the data points
        plot((eval(c4)),(eval(c5)), '.', 'Color', color)
        hold on

        % Plot the regression line
        %h = plot(mdl);
        plot(x_fit, y_fit, 'k', 'LineWidth', 1);
       % h(2).Color = 'k';  % Set regression line color
        %h(2).LineWidth = 1;
       % h(1).Visible = 'off';  % Hide original data points in fitlm plot
      %  h(3).Visible = 'off';  % Hide confidence bounds
      %  h(4).Visible = 'off';  % Hide confidence bounds

        % Customize the plot
        %legend(h, 'off')
        %legend boxoff
        %delete(h(1))%marker
       % delete(h(3))%bound
        %delete(h(4))%bound
        xlabel('');  % Remove x-axis label
        ylabel('');  % Remove y-axis label
        title('');

        % Add text with correlation and p-value
        text(-750, 700, ['PRCC =  ' num2str(round(r, 2)) ' ,  P-value \approx '  num2str(round(p, 2))], 'fontsize', 7)
        
        % Set axis limits and ticks
        xlim([-800, 800])
        ylim([- 800, 800])
        xticks([-1000, -500, 0, 500, 1000])
        yticks([-1000, -500, 0, 500, 1000])
        
        % Set figure properties for saving
        set(gca, 'FontSize', 8)
        set(gcf, 'PaperPositionMode', 'manual'); 
        set(gcf, 'PaperUnits', 'centimeters'); 
        set(gcf, 'PaperPosition', [0 0  4.5 3.5]); 
        set(gcf, 'PaperSize', [4.5 3.5]);

    % Ensure 'metric' is a character vector
    if isstring(metric)
        metric = char(metric);
    elseif ~ischar(metric)
        metric = num2str(metric);
    end

        % Save the figure as a PDF
        pathname = fileparts('/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new/Plot/');
        figfile = fullfile(pathname, [num2str(i) '_' metric '.png']);
        % Save the figure as PNG with high resolution (1080 DPI)
        print(fig(i), figfile, '-dpng', '-r1080');
        hold off
    
        hold off  % Release hold after plotting
    end
end
end
