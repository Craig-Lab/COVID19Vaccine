close all
clc
clear variables

file_to_load = 0;
%0 - 3000 data point file on Aug 21

alpha = 0.01;

switch file_to_load
    case 0
        file_title='prcc_Model_LHS';
       % load([file_title(6:end),'.mat']);
        load('Model_LHS.mat');
        
end

%L
[prcc{1},sign{1},sign_label{1}]=PRCC(LHSmatrix,peak_A1,1,PRCC_var,alpha);
%V
[prcc{2},sign{2},sign_label{2}]=PRCC(LHSmatrix,peak_A2,1,PRCC_var,alpha);
%Th
[prcc{3},sign{3},sign_label{3}]=PRCC(LHSmatrix,peak_time1,1,PRCC_var,alpha);
%B
[prcc{4},sign{4},sign_label{4}]=PRCC(LHSmatrix,peak_time2,1,PRCC_var,alpha);
%GB
[prcc{5},sign{5},sign_label{5}]=PRCC(LHSmatrix,auc_A1,1,PRCC_var,alpha);
%SP
[prcc{6},sign{6},sign_label{6}]=PRCC(LHSmatrix,auc_A2,1,PRCC_var,alpha);
%LP
[prcc{7},sign{7},sign_label{7}]=PRCC(LHSmatrix,auc_total,1,PRCC_var,alpha);
%M
%[prcc{8},sign{8},sign_label{8}]=PRCC(LHSmatrix,R,1:length(time_points),PRCC_var,alpha);
%A
%[prcc{9},sign{9},sign_label{9}]=PRCC(LHSmatrix,R,1:length(time_points),PRCC_var,alpha);
%I
%[prcc{10},sign{10},sign_label{10}]=PRCC(LHSmatrix,R,1:length(time_points),PRCC_var,alpha);

save(file_title,'prcc','sign','sign_label','PRCC_var');



