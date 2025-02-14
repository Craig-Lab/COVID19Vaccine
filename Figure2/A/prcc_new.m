function PRCCs = prcc_new(varargin)

%% Description of Function, Inputs, and Outputs
% This function serves as a means of calculating the partial rank
% correlation coefficients (PRCCs) of a set of equally sized (n,1) input 
% arrays, p1, p2,... where the last input is the response of the model for 
% each iteration of the given parameters. The output of this function is a 
% vector of the PRCC for each parameter provided, with a corresponding bar
% plot. For optimal usage, ensure an appropriately large number of
% iterations for each parameter and model response. Stratified sampling 
% techniques are suggested for obtaining these iterations of parameters, 
% so as to limit the computational burden of this function.

%% PRCC and Applications
% Most practical applications of partial rank correlation coefficient
% analysis take place in uncertainty and sensitivity analysis,
% particularly in biological systems with comparatively high degrees of
% parameter uncertainty. These are used as a means of exploring a system's 
% parameter space while attempting to control for the impacts of other
% parameters on the response, allowing each of them to be analyzed 
% simultaneously. This is typically used, in comparison to a partial 
% correlation coefficient test, in circumstances of nonlinear but monotonic 
% relationships between a function and its parameters. For further
% documentation, see:
%
% Marino, Simeone, et al. "A methodology for performing global uncertainty 
% and sensitivity analysis in systems biology." Journal of theoretical
% biology 254.1 (2008): 178-196.
%
% Gomero, Boloye. “Latin Hypercube Sampling and Partial Rank Correlation 
% Coefficient Analysis Applied to an Optimal Control Problem.” (2012).

%% Development and Author Information
%
% * Developed in Matlab 2022a, expected compatability with all versions
%   later than Matlab 2006a, no additional toolboxes required for usage
%
% * Tristen M Jackson, Florida State University, tjackson@math.fsu.edu
%       All rights reserved
%
% * Last Modified: October 14, 2022

i=1;
j=1;
k=[1:1:100000];


    if nargin == 0
        error('No input arguments')
    end

    if nargin == 1
        error('Not enough input arguments. Inputs must include at least two vectors')
    end

    for i = 2:nargin
        if length(varargin{i})==length(varargin{i-1})
        else error('all input arrays must be the same length')
        end
    end


    if nargin > 1
        for j=1:nargin
            Filler{:,j} = sort([varargin{j}], 'descend');
        end   

        Parameters=Filler(:,1:end-1);
        END = Filler(:,end);
        EachColumnRemoved = [];

        for col = 1:size(Parameters,2)
            EachColumnRemoved{end+1} = [Parameters(:,1:col-1), Parameters(:,col+1:end)];
        end

        for i=1:nargin-1
            [~,~, ResidualsEachColumnRemoved{1,i}] = regress(END{1,1}, cell2mat(EachColumnRemoved{1,i}));
            [~,~, ResidualsEachColumn{1,i}] = regress(cell2mat(Parameters(1,i)), cell2mat(EachColumnRemoved{1,i}));
            Coefficients{1,i} = corrcoef(ResidualsEachColumn{1,i}, ResidualsEachColumnRemoved{1,i});
            [PRCCs(1,i)] = Coefficients{1,i}(1,2);
        end

    bar(diag(PRCCs), 'stacked')
    legendNames = "p" + k;
    legend(legendNames);
    xlabel('Parameters')
    title('Partial Rank Correlation Coefficients')
    end
end