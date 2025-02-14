% lipid nanoparticles
p.d_l = 0.1;
p.del_lv = 0.15; %0.246 -0.1

% vaccinated cells
p.d_v = 0.89;%fitted


% T follicular helper cells
p.del_tv = 0.73;%4.98;%0.11; %4.98;%0.035;%4.98;%%%%%%%%%%%%4.98
p.d_t = 0.57;%fitted

    
% activated b cells
p.del_bt = 0.36;%fitted
p.rho_g =0.69 ;%generation rate of GC B cells
p.rho_s =0.6 ;% generation rate of SLPCs ( p.del_sp)
p.d_b = 0.02;%half-life 6 days

% germinal center B cells
p.del_ig = 0.72; %binding rate of IL-21
p.beta_g = 2.77; % fitted 0.057;% 0.0177;%0.079; %8;%0.26;%0.01;%1; %0.0046*1e3; %9/24; %per cell per day
p.p_g = 0.1;
p.d_g = 0.49;%fitted 

% long-lived plasma cells
p.p_p = 0.27;%0.9;
p.d_p =  0.0112;% mean life 180 days 0.0027;%0.042;% 360 days 38; %half life 18 days 0.042;%71;%half life 11.6 chapin's paper %0.071; 

% short-lived plasma cells
p.d_s = 0.29; %lifespan 5 days fitted
    
% memory B cells
p.d_m = 0.0029;
p.p_m = 0.85; 
p.beta_m = 0.009;%lower bound of literature
%p.SV = 0.059;%
p.ht = 1;%%%%%
p.p_p2 = 0.9;
% antibody
p.alpha_p = 18; %0.0424antibody production rate by LLPCs
p.alpha_s =272;%0.0086antibody production rate by SLPCs
p.d_a =  0.033; %halflife 60 days; mean 90

% interleukin
p.rho_i =2.07;
p.d_i = 3.47;%;52.63;%
p.SI = 0.03;%fitted


p.L0 =  1;%1; %30 ug/30mL pfizer 100ug moderna
p.V0 = 0;
p.Th0=0;
p.B0 = 0;%avtivated b cells
p.GB0 = 0;
p.LP0 = 0;
p.SP0 = 0;
p.M0=0;
p.A0=0;%22.5;
p.I0=0;%5.18