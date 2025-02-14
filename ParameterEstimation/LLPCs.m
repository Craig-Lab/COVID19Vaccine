function [time,sol,solstruc] = LLPCs(p,tspan)
p.IC = [1;0;0;0;0;0;0;0];
%p.IC = [p.L0;p.V0;p.Th0;p.B0;p.GB0;p.SP0;p.I0];

opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
%solstruc = ddesd_f5(@ddefun,@(t,y) delayP(t,y,p),@history,tspan,opts);
solstruc = ode45(@odefun,tspan,p.IC,opts);

time = solstruc.x;%linspace(tspan(1),tspan(end),1000);
sol = solstruc.y;%deval(solstruc,time);

%------------------------------------------------------------------------


function dydt = odefun(t,y)
   
    L = y(1);%lipid nanoparticles
    V = y(2);%vaccinated cells
    Th = y(3);%helper T cells
    B = y(4);
    GB = y(5);%Germinal center GB cells
    SP = y(6);
    LP = y(7);
    I= y(8); %interleukin

  
    dL = -p.del_lv *L- p.d_l *L;
    dV = p.del_lv *L-p.d_v*V ;
    dTh =p.del_tv*V-p.d_t*Th;
    dB = p.del_bt*Th - p.d_b*B-B*p.sum;%(p.rho_g+p.rho_s);
    dGB = B*p.rho_g+p.beta_g*(2*p.p_g-1)*p.del_ig*I/(p.SI+I)*GB-p.d_g*GB;
    dSP = B*p.rho_s-p.d_s*SP;
    dLP = (1-p.p_g)*p.p_p*p.beta_g*GB-p.d_p*LP;
    dI = p.rho_i*Th - p.d_i*I;
    
    dydt = [dL;dV;dTh;dB;dGB;dSP;dLP;dI];

end

end


