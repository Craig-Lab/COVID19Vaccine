function [time,sol,solstruc] = ActivatedBcells(p,tspan)

%p.IC = [p.L0;p.V0;p.Th0;p.B0];
p.IC = [1e3;0;0;0];
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
solstruc = ode45(@odefun,tspan,p.IC,opts);

time = solstruc.x;%linspace(tspan(1),tspan(end),1000);
sol = solstruc.y;%deval(solstruc,time);

%------------------------------------------------------------------------


function dydt = odefun(t,y)
   
    L = y(1);%lipid nanoparticles
    V = y(2);%vaccinated cells
    Th = y(3);%helper T cells
    B = y(4);
   
    dL = -p.del_lv *L- p.d_l *L;
    dV = p.del_lv *L-p.d_v*V ;
    dTh =p.del_tv*V-p.d_t*Th;
    dB = p.del_bt*Th - p.d_b*B-B*p.a;
 
    dydt = [dL;dV;dTh;dB];

end

end
