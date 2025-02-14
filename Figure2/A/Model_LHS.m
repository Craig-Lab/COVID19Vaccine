%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below
clear variables;
close all;

if exist('Model_LHS_10000_Aug2120','file')~=0
    resp = input('This file already exists, rerun? Y/N\n','s');
    if strcmpi(resp,'Y')
        run_sims=1;
    else
        run_sims=0;
    end
else
    run_sims = 1;
end

if run_sims

    %% Sample size N
    runs=1000;

    %% LHS MATRIX  %%
    Parameter_settings_LHS;

    fprintf('Setting random parameters...\n');

    del_lv_LHS=LHS_Call(0.5*params.del_lv, params.del_lv, 1.5*params.del_lv, 0, runs,'unif');%1
    d_l_LHS=LHS_Call(0.5*params.d_l, params.d_l, 1.5*params.d_l, 0, runs,'unif');%2
    d_v_LHS=LHS_Call(0.5*params.d_v, params.d_v, 1.5*params.d_v, 0, runs,'unif');%3
    del_tv_LHS=LHS_Call(0.5*params.del_tv, params.del_tv, 1.5*params.del_tv, 0, runs,'unif');%4
    d_t_LHS=LHS_Call(0.5*params.d_t , params.d_t , 1.5*params.d_t , 0, runs,'unif');%5
    del_bt_LHS=LHS_Call(0.5*params.del_bt , params.del_bt , 1.5*params.del_bt , 0, runs,'unif');%6
    rho_g_LHS=LHS_Call(0.5*params.rho_g , params.rho_g , 1.5*params.rho_g , 0, runs,'unif');%7
    rho_s_LHS=LHS_Call(0.5*params.rho_s , params.rho_s , 1.5*params.rho_s , 0, runs,'unif');%8
    d_b_LHS=LHS_Call(0.5*params.d_b , params.d_b , 1.5*params.d_b , 0, runs,'unif');%9 
    del_ig_LHS=LHS_Call(0.5*params.del_ig , params.del_ig , 1.5*params.del_ig , 0, runs,'unif');%10
    beta_g_LHS=LHS_Call(0.5*params.beta_g , params.beta_g , 1.5*params.beta_g , 0, runs,'unif');%11
    p_g_LHS=LHS_Call(0.5*params.p_g  , params.p_g  , 1.5*params.p_g  , 0, runs,'unif');%12
    d_g_LHS=LHS_Call(0.5*params.d_g  , params.d_g , 1.5*params.d_g  , 0, runs,'unif');%13
    p_p_LHS=LHS_Call(0.5*params.p_p  , params.p_p  , 1.5*params.p_p  , 0, runs,'unif');%14
    d_p_LHS=LHS_Call(0.5*params.d_p  , params.d_p  , 1.5*params.d_p , 0, runs,'unif');%15
    d_s_LHS=LHS_Call(0.5*params.d_s  , params.d_s , 1.5*params.d_s  , 0, runs,'unif');%16
    d_m_LHS=LHS_Call(0.5*params.d_m  , params.d_m  , 1.5*params.d_m  , 0, runs,'unif');%17
    p_m_LHS=LHS_Call(0.5*params.p_m  , params.p_m  , 1.5*params.p_m  , 0, runs,'unif');%18
    beta_m_LHS=LHS_Call(0.5*params.beta_m  , params.beta_m  , 1.5*params.beta_m  , 0, runs,'unif');%19
    p_p2_LHS=LHS_Call(0.5*params.p_p2  , params.p_p2   , 1.5*params.p_p2,  0, runs,'unif');%20
    ht_LHS=LHS_Call(0.5*params.ht  , params.ht  , 1.5*params.ht  , 0, runs,'unif');%21
    alpha_p_LHS=LHS_Call(0.5*params.alpha_p  , params.alpha_p  , 1.5*params.alpha_p  , 0, runs,'unif');%22
    alpha_s_LHS=LHS_Call(0.5*params.alpha_s  , params.alpha_s  , 1.5*params.alpha_s  , 0, runs,'unif');%23
    d_a_LHS=LHS_Call(0.5*params.d_a  , params.d_a  , 1.5*params.d_a  , 0, runs,'unif');%24
    rho_i_LHS=LHS_Call(0.5*params.rho_i  , params.rho_i  , 1.5*params.rho_i  , 0, runs,'unif');%25
    d_i_LHS=LHS_Call(0.5*params.d_i  , params.d_i  , 1.5*params.d_i  , 0, runs,'unif');%26
    SI_LHS=LHS_Call(0.5*params.SI  , params.SI  , 1.5*params.SI  , 0, runs,'unif');%27
     
    fprintf('Random parameters set...\n');

    %% LHS MATRIX and PARAMETER LABELS
        LHSmatrix=[del_lv_LHS,d_l_LHS,d_v_LHS,del_tv_LHS,d_t_LHS,del_bt_LHS,rho_g_LHS,rho_s_LHS,d_b_LHS,...
            del_ig_LHS, beta_g_LHS,p_g_LHS,d_g_LHS,p_p_LHS,d_p_LHS,d_s_LHS,d_m_LHS,p_m_LHS,beta_m_LHS,...
            p_p2_LHS,ht_LHS,alpha_p_LHS,alpha_s_LHS,d_a_LHS,rho_i_LHS,d_i_LHS,SI_LHS];
    tic
        peak_A1 = zeros(1,runs);
        peak_A2 = zeros(1,runs);
        peak_time1 =zeros(1,runs);
        peak_time2 =zeros(1,runs);
        auc_A1 = zeros(1,runs);
        auc_A2 = zeros(1,runs);
        auc_total =zeros(1,runs);

    for x=1:runs %Run solution x times choosing different values
       
        %f=@ODE_LHS;
        fprintf('Running simulation for parameter set %i...\n',x);
        params.L0 = 1;
        params.V0 = 0;
        params.Th0 = 0;
        params.B0 =  0;
        params.GB0 =  0;
        params.LP0 =  0;
        params.SP0 = 0;
        params.M0= 0;
        params.A0= 0;
        params.I0 =  0;
        params.del_lv = LHSmatrix(x,1);
        params.d_l = LHSmatrix(x,2);
        params.d_v = LHSmatrix(x,3);
        params.del_tv = LHSmatrix(x,4);
        params.d_t = LHSmatrix(x,5);
        params.del_bt = LHSmatrix(x,6);
        params.rho_g = LHSmatrix(x,7);
        params.rho_s = LHSmatrix(x,8);
        params.d_b = LHSmatrix(x,9);
        params.del_ig= LHSmatrix(x,10);
        params.beta_g = LHSmatrix(x,11);
        params.p_g = LHSmatrix(x,12);
        params.d_g =  LHSmatrix(x,13);
        params.p_p= LHSmatrix(x,14);
        params.d_p = LHSmatrix(x,15);
        params.d_s = LHSmatrix(x,16);
        params.d_m = LHSmatrix(x,17);
        params.p_m = LHSmatrix(x,18);
        params.beta_m = LHSmatrix(x,19);
        params.p_p2 = LHSmatrix(x,20);
        params.ht = LHSmatrix(x,21);
        params.alpha_p = LHSmatrix(x,22);
        params.alpha_s = LHSmatrix(x,23);
        params.d_a =  LHSmatrix(x,24);
        params.rho_i = LHSmatrix(x,25);
        params.d_i = LHSmatrix(x,26);
        params.SI = LHSmatrix(x,27);
        [time1,sol1] = Humoral_response_model1(params,tspan);
        
        time1_deval= linspace(0,60,1e3);
        sol1_deval = interp1(time1,sol1',time1_deval);
        A1=[sol1_deval,time1_deval']; 
        params.L0 = sol1_deval(1000,1)+1;
        params.V0 = sol1_deval(1000,2);
        params.Th0 = sol1_deval(1000,3);
        params.B0 =  sol1_deval(1000,4);
        params.GB0 =  sol1_deval(1000,5);
        params.LP0 =  sol1_deval(1000,6);
        params.SP0 = sol1_deval(1000,7);
        params.M0= sol1_deval(1000,8);
        params.A0= sol1_deval(1000,9);
        params.I0 =  sol1_deval(1000,10);
        params.IC = [params.L0;params.V0;params.Th0;params.B0;params.GB0;...
            params.LP0;params.SP0;params.M0;params.A0;params.I0];
        tspan2 = [60, 270];
      
       params.lam2=0;
        [time2,sol2] = Humoral_response_model2(params,tspan2);
        time2_deval= linspace(60,270,1e3);
        sol2_deval = interp1(time2,sol2',time2_deval);
        A2=[sol2_deval,time2_deval' ]; % [time y]
        A = [A1;A2];
      
        [M1, I1] = max(A1(:,9));
        [M2, I2] = max(A2(:,9));
        peak_A1(x) = M1;
        peak_A2(x) = M2;
        peak_time1(x) = time1_deval(I1);
        peak_time2(x) = time2_deval(I2);
        auc_A1(x) = trapz(time1,sol1(9,:)*1e3);
        auc_A2(x) = trapz(time2,sol2(9,:)*1e3);
        auc_total(x) = trapz(time1,sol1(9,:)*1e3)+ trapz(time2,sol2(9,:)*1e3);
 
     

        %% Save only the outputs at the time points of interest [time_points]:
        %% MORE EFFICIENT
        L_lhs(:,x)=A(time_points,1);
        V_lhs(:,x)=A(time_points,2);
        Th_lhs(:,x)=A(time_points,3);
        B_lhs(:,x)=A(time_points,4);
        GB_lhs(:,x)=A(time_points,5);
        LP_lhs(:,x)=A(time_points,6);
        SP_lhs(:,x)=A(time_points,7);
        M_lhs(:,x)=A(time_points,8);
        A_lhs(:,x)=A(time_points,9);
        I_lhs(:,x)=A(time_points,10);

    end
    time_end = toc;
    fprintf('Simulations completed in %3.2fs...\n',time_end);
    %% Save the workspace
    save Model_LHS.mat;
end

