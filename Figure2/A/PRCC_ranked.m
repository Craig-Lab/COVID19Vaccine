%% It calculates PRCCs and their significances
%% (uncorrected p-value, Bonferroni correction and 
%% Benjamini and Hochberg False Discovery Rate correction)
%% LHSmatrix: LHS matrix (N x k) %%
%% Y: output matrix (time x N) %%
%% s: time points to test (row vector) %%
%% PRCC_var: vector of strings {'p1','p2',...,'pk'}
%%            to label all the parameters varied in the LHS
%% by Simeone Marino, May 29 2007 %%
%% N.B.: this version uses ONE output at a time

function [prcc,sign,sign_label]=PRCC_ranked(X,Y,s,PRCC_var)

%LHSmatrix=LHS; % Define LHS matrix
Y=Y(s,:)';% Define the output. Comment out if the Y is already 
          % a subset of all the time points and it already comprises
          % ONLY the s rows of interest
[a,k]=size(X); % Define the size of LHS matrix
[b,out]=size(Y);
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
   [rho,p]=corr(eval(c4),eval(c5));
    for j=1:out
        c6=['prcc_',num2str(i),'(',num2str(j),')=rho(1,',num2str(j),');'];
        c7=['prcc_sign_',num2str(i),'(',num2str(j),')=p(1,',num2str(j),');'];
        eval(c6);
        eval(c7);
    end
    c8=['clear r',num2str(i),';'];
    eval(c8);
end
prcc=[];
prcc_sign=[];
for i=1:k
    d1=['prcc=[prcc ; prcc_',num2str(i),'];'];
    eval(d1);
    d2=['prcc_sign=[prcc_sign ; prcc_sign_',num2str(i),'];'];
    eval(d2);
end
[length(s) k out];
PRCCs=prcc';
uncorrected_sign=prcc_sign';
prcc=PRCCs;
sign=uncorrected_sign;

%% Multiple tests correction: Bonferroni and FDR
%tests=length(s)*k; % # of tests performed
%correction_factor=tests;
%Bonf_sign=uncorrected_sign*tests;
%sign_new=[];
%for i=1:length(s)
%    sign_new=[sign_new;(1:k)',ones(k,1)*s(i),sign(i,:)'];
%end
%sign_new=sortrows(sign_new,3);
%for j=2:k
%    sign_new(j,3)=sign_new(j,3)*(tests/(tests-j+1));
%end
%sign_new=sortrows(sign_new,[2 1]); % FDR correction
%sign_new=sign_new(:,3)';
%for i=1:length(s)
%    FDRsign(i,:)=[sign_new(1+k*(i-1):i*k)];
%end
%uncorrected_sign; % uncorrected p-value
%Bonf_sign;  % Bonferroni correction
%FDRsign; % FDR correction

sign_label_struct=struct;
sign_label_struct.uncorrected_sign=uncorrected_sign;
sign_label=sign_label_struct;