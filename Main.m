clear all; close all; clc;

constants

%% ---> Time discretization
n_periodes = 1;
Time = 1000;
dt = 0.02;
dtt = 0.02;

tEnd    = Time*n_periodes;     
nStep   = tEnd/dt;      %---> Number of global time steps
nnStep  = tEnd/dtt;     %---> Number of intern time steps

%% ---> Initial conditions
V0 = -86;
c_na = 11.6; c_k = 138.3; c_ca = 0.08*10^(-3); c_srca = 0.56;
g_m = 0; g_h = 0.75; g_j = 0.75; g_d = 0; g_f = 1; g_fca = 1; g_r = 0; g_s = 1; g_xs = 0; g_xr1 = 0; g_xr2 = 0; g_k1 = 0.5; g_g = 1;

%% ---> Iterative process
% Potential vector
V = zeros(nStep,1); V(1,1)= V0; V_ant = V(1,1);

% Gates vectors
G_m = zeros(nStep,1); G_m(1,1) = g_m; gant_m = g_m;
G_h = zeros(nStep,1); G_h(1,1) = g_h; gant_h = g_h;
G_j = zeros(nStep,1); G_j(1,1) = g_j; gant_j = g_j;
G_d = zeros(nStep,1); G_d(1,1) = g_d; gant_d = g_d;
G_f = zeros(nStep,1); G_f(1,1) = g_f; gant_f = g_f;
G_fca = zeros(nStep,1); G_fca(1,1) = g_fca; gant_fca = g_fca;
G_r = zeros(nStep,1); G_r(1,1) = g_r; gant_r = g_r;
G_s = zeros(nStep,1); G_s(1,1) = g_s; gant_s = g_s;
G_xs = zeros(nStep,1); G_xs(1,1) = g_xs; gant_xs = g_xs;
G_xr1 = zeros(nStep,1); G_xr1(1,1) = g_xr1; gant_xr1 = g_xr1;
G_xr2 = zeros(nStep,1); G_xr2(1,1) = g_xr2; gant_xr2 = g_xr2;
G_k1 = zeros(nStep,1); G_k1(1,1) = g_k1; gant_k1 = g_k1;
G_g = zeros(nStep,1); G_g(1,1) = g_g; gant_g = g_g;

% Concentration vectors
C_k = zeros(nStep,1); C_k(1,1) = c_k; cant_k = c_k;
C_na = zeros(nStep,1); C_na(1,1) = c_na; cant_na = c_na;
C_ca = zeros(nStep,1); C_ca(1,1) = c_ca; cant_ca = c_ca;
C_srca = zeros(nStep,1); C_srca(1,1) = c_srca; cant_srca = c_srca;

for n = 2:nStep + 1
    n
    
    % Update gateI gate variables (10)
    gnou_m = gateI_m(V_ant,gant_m,dt); G_m(n,1) = gnou_m;
    gnou_h = gateI_h(V_ant,gant_h,dt); G_h(n,1) = gnou_h;
    gnou_j = gateI_j(V_ant,gant_j,dt); G_j(n,1) = gnou_j;
    gnou_xr1 = gateI_xr1(V_ant,gant_xr1,dt); G_xr1(n,1) = gnou_xr1;
    gnou_xr2 = gateI_xr2(V_ant,gant_xr2,dt); G_xr2(n,1) = gnou_xr2;
    gnou_xs = gateI_xs(V_ant,gant_xs,dt); G_xs(n,1) = gnou_xs;
    gnou_r = gateI_r(V_ant,gant_r,dt); G_r(n,1) = gnou_r;
    gnou_s = gateI_s(V_ant,gant_s,dt); G_s(n,1) = gnou_s;
    gnou_d = gateI_d(V_ant,gant_d,dt); G_d(n,1) = gnou_d;
    gnou_f = gateI_f(V_ant,gant_f,dt); G_f(n,1) = gnou_f;
    
    % Iterations for c_ion, gateII and I_ion variables
    [vect_x,vect_r,iter_fi,gnou_k1,gnou_fca,gnou_g,JRES0,RES0] = newton_iter_local(nnStep,10^(-9),cant_k,cant_na,cant_ca,cant_srca,...
            gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,gnou_r,gnou_s,gnou_d,gnou_f,gant_fca,gant_g,V_ant,dtt);
    
    cnou_k = vect_x(1,iter_fi); C_k(n,1) = cnou_k;
    cnou_na = vect_x(2,iter_fi); C_na(n,1) = cnou_na;
    cnou_ca = vect_x(3,iter_fi); C_ca(n,1) = cnou_ca;
    cnou_srca = vect_x(4,iter_fi); C_srca(n,1) = cnou_srca;
    G_k1(n,1) = gnou_k1;
    G_fca(n,1) = gnou_fca;
    G_g(n,1) = gnou_g;
      
    % Global time iterations
    V_nou = newton_iter_time(V_ant,gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,gnou_r,gnou_s,gnou_d,gnou_f,gnou_k1,gnou_fca,cnou_na,cnou_k,cnou_ca,dt,n);
    V(n,1) = V_nou;
    
    % Update old (ant) variables
    V_ant = V_nou;
    cant_k = cnou_k;
    cant_na = cnou_na;
    cant_ca = cnou_ca;
    cant_srca = cnou_srca;
    gant_m = gnou_m;
    gant_h = gnou_h;
    gant_j = gnou_j;
    gant_xr1 = gnou_xr1;
    gant_xr2 = gnou_xr2;
    gant_xs = gnou_xs;
    gant_r = gnou_r;
    gant_s = gnou_s;
    gant_d = gnou_d;
    gant_f = gnou_f;
    gant_k1 = gnou_k1;
    gant_fca = gnou_fca;
    gant_g = gnou_g;
        
end

%% ---> Plot: membrane potential and concentrations
figure(1)
t_vect = 0:dt:tEnd;
y1 = V(:,1);

plot(t_vect.*0.001,y1,'b')
title('Membrane potential')
xlabel('s') 
ylabel('mV')

%% ---> Plot: Concentrations
figure(2)
y2 = C_k(:,1);
y3 = C_na(:,1);
y4 = C_ca(:,1);
y5 = C_srca(:,1);

ax2 = subplot(2,2,1);
plot(ax2,t_vect.*0.001,y2,'r')
title('c K')
xlabel('s') 
ylabel('mM')

ax3 = subplot(2,2,2);
plot(ax3,t_vect.*0.001,y3,'r')
title('c Na')
xlabel('s') 
ylabel('mM')

ax4 = subplot(2,2,3);
plot(ax4,t_vect.*0.001,y4*1000,'r')
title('c Ca')
xlabel('s') 
ylabel('microM')

ax5 = subplot(2,2,4);
plot(ax5,t_vect.*0.001,y5,'r')
title('c Casr')
xlabel('s') 
ylabel('mM')


%% ---> Plot: voltage dependent evolution of tau and ginf
figure(3)
x = V(:,1);

% gateI_m
ginf_m = (1+exp((-56.86-x)./9.03)).^(-2);
tau_m = 0.1.*(1+exp((-60-x)./5)).^(-1).*((1+exp((x+35)/5)).^(-1)+(1+exp((x-50)./200)).^(-1));
plot(subplot(5,4,1),x,tau_m,'c')
title('tau m')
xlabel('mV')
plot(subplot(5,4,2),x,ginf_m,'m')
title('ginf m')
xlabel('mV')

% gateI_h gateI_j
ginf_h = (1+exp((x+71.55)./7.43)).^(-2);
tau_h = 0.1688.*(1+exp(-(V+10.66)./11.1)).*(V >= -40) + (0.057.*exp(-(80+V)./6.8)+2.7.*exp(0.079.*V)+3.1.*(10^5).*exp(0.3485.*V)).^(-1).*(V < -40);
plot(subplot(5,4,3),x,tau_h,'c')
title('tau h')
xlabel('mV')
plot(subplot(5,4,4),x,ginf_h,'m')
title('ginf h i ginf j')
xlabel('mV')
tau_j = (0.6.*exp(0.057.*V).*(1+exp(-0.1.*(V+32))).^(-1)).^(-1).*(V >= -40)...
            + ((-2.5428.*(10^4).*exp(0.2444.*V)-6.948.*(10^(-6)).*exp(-0.04391.*V)).*(V+37.78).*(1+exp(0.311.*(V+79.23))).^(-1)+0.02424.*exp(-0.01052.*V).*(1+exp(-0.1378.*(V+40.14))).^(-1)).^(-1).*(V < -40);
plot(subplot(5,4,5),x,tau_j,'c')
title('tau j')
xlabel('mV')

% gateII_fca
ginf_fca = 0.685.*((1+(C_ca./0.000325).^8).^(-1)+0.1.*(1+exp((C_ca-0.0005)./0.0001)).^(-1)+0.2.*((1+exp((C_ca-0.00075)./0.0008))).^(-1) +0.23);
plot(subplot(5,4,6),C_ca,ginf_fca,'m')
title('ginf fca')
xlabel('c Ca [microM]')

% gateI_f
ginf_f = (1+exp((x+20)./7)).^(-1);
tau_f = (1125.*exp(-(x+27).^2./240))+(165.*(1+exp((25-x)./10)).^(-1))+80;
plot(subplot(5,4,7),x,tau_f,'c')
title('tau f')
xlabel('mV')
plot(subplot(5,4,8),x,ginf_f,'m')
title('ginf f')
xlabel('mV')

% gateI_d
ginf_d = (1+exp((-5-x)./7.5)).^(-1);
tau_d = ((1.4).*(1+exp((-35-x)./13)).^(-1)+0.25).*((1.4).*(1+exp((x+5)./5)).^(-1))+(1+exp((50-x)./20)).^(-1);
plot(subplot(5,4,9),x,tau_d,'c')
title('tau d')
xlabel('mV')
plot(subplot(5,4,10),x,ginf_d,'m')
title('ginf d')
xlabel('mV')

% gateI_r
ginf_r = (1+exp((20-x)./6)).^(-1);
tau_r = 9.5.*exp(-(x+40).^2./1800)+0.8;
plot(subplot(5,4,11),x,tau_r,'c')
title('tau r')
xlabel('mV')
plot(subplot(5,4,12),x,ginf_r,'m')
title('ginf r')
xlabel('mV')

% gateI_s
ginf_s = (1+exp((x+20)./5)).^(-1);
tau_s = 85.*exp(-(x+45).^2./320)+5.*(1+exp((x-20)./5)).^(-1)+3;
plot(subplot(5,4,13),x,tau_s,'c')
title('tau s')
xlabel('mV')
plot(subplot(5,4,14),x,ginf_s,'m')
title('ginf s')
xlabel('mV')

% gateI_xs
ginf_xs = (1+exp((-5-x)./14)).^(-1);
tau_xs = 1100.*(1+exp((-10-x)./6)).^(-1/2).*(1+exp((x-60)./20)).^(-1);
plot(subplot(5,4,15),x,tau_xs,'c')
title('tau xs')
xlabel('mV')
plot(subplot(5,4,16),x,ginf_xs,'m')
title('ginf xs')
xlabel('mV')

% gateI_xr1
ginf_xr1 = (1+exp((-26-x)./7)).^(-1);
tau_xr1 = 2700.*(1+exp((-45-x)./10)).^(-1).*((1+exp((x+30)./11.5))).^(-1);
plot(subplot(5,4,17),x,tau_xr1,'c')
title('tau xr1')
xlabel('mV')
plot(subplot(5,4,18),x,ginf_xr1,'m')
title('ginf xr1')
xlabel('mV')

% gateI_xr2
ginf_xr2 = (1+exp((x+88)./24)).^(-1);
tau_xr2 = 3.36.*(1+exp((-60-x)./20)).^(-1).*(1+exp((x-60)./20)).^(-1);
plot(subplot(5,4,19),x,tau_xr2,'c')
title('tau xr2')
xlabel('mV')
plot(subplot(5,4,20),x,ginf_xr2,'m')
title('ginf xr2')
xlabel('mV')


%% ---> Plot: Gates
figure(4)
plot(subplot(4,4,1),t_vect.*0.001,G_m)
title('g m')
xlabel('s')
plot(subplot(4,4,2),t_vect.*0.001,G_h)
title('g h')
xlabel('s')
plot(subplot(4,4,3),t_vect.*0.001,G_j)
title('g j')
xlabel('s')
plot(subplot(4,4,4),t_vect.*0.001,G_d)
title('g d')
xlabel('s')
plot(subplot(4,4,5),t_vect.*0.001,G_f)
title('g f')
xlabel('s')
plot(subplot(4,4,6),t_vect.*0.001,G_fca)
title('g fca')
xlabel('s')
plot(subplot(4,4,7),t_vect.*0.001,G_r)
title('g r')
xlabel('s')
plot(subplot(4,4,8),t_vect.*0.001,G_s)
title('g s')
xlabel('s')
plot(subplot(4,4,9),t_vect.*0.001,G_xs)
title('g xs')
xlabel('s')
plot(subplot(4,4,10),t_vect.*0.001,G_xr1)
title('g xr1')
xlabel('s')
plot(subplot(4,4,11),t_vect.*0.001,G_xr2)
title('g xr2')
xlabel('s')
plot(subplot(4,4,12),t_vect.*0.001,G_k1)
title('g k1')
xlabel('s')
plot(subplot(4,4,13),t_vect.*0.001,G_g)
title('g g')
xlabel('s')

%% ---> Plot: Currents
figure(5)
phi_na = R*T*log(c_na0./C_na)/(z_na*F);
phi_k = R*T*log(c_k0./C_k)/(z_k*F);
phi_ca = R*T*log(c_ca0./C_ca)/(z_ca*F);
phi_ks = R*T/F*log((c_k0+p_kna*c_na0).*(C_k+p_kna.*C_na).^(-1));

I_na = Cmax_na.*G_m.^3.*G_h.*G_j.*(V-phi_na);
plot(subplot(4,4,1),t_vect.*0.001,I_na)
title('I na')
xlabel('s') 
ylabel('pA/pF') 

I_bna = Cmax_bna.*(V-phi_na);
plot(subplot(4,4,2),t_vect.*0.001,I_bna)
title('I bna')
xlabel('s') 
ylabel('pA/pF')

I_nak=Imax_nak*(c_k0.*C_na).*((C_na+c_nak).*(c_k0+c_kna).*(1+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T)))).^(-1);
plot(subplot(4,4,3),t_vect.*0.001,I_nak)
title('I nak')
xlabel('s') 
ylabel('pA/pF')

I_naca=Imax_naca*(exp(y.*V*F/(R*T)).*C_na.^3*c_ca0-exp((y-1).*V*F/(R*T))*c_na0^3.*C_ca*y_naca)...
            .*((c_naca^3+c_na0^3)*(c_cana+c_ca0)*(1+k_naca*exp((y-1).*V*F/(R*T)))).^(-1);
plot(subplot(4,4,4),t_vect.*0.001,I_naca)
title('I naca')
xlabel('s') 
ylabel('pA/pF')

I_k1 = Cmax_k1*G_k1.*(c_k0./5.4).^(1/2).*(V-phi_k);
plot(subplot(4,4,5),t_vect.*0.001,I_k1)
title('I k1')
xlabel('s') 
ylabel('pA/pF')

I_kr=Cmax_kr.*G_xr1.*G_xr2.*(c_k0/(5.4))^(1/2).*(V-phi_k);
plot(subplot(4,4,6),t_vect.*0.001,I_kr)
title('I kr')
xlabel('s') 
ylabel('pA/pF')

I_ks = Cmax_ksepi.*G_xs.^2.*(V-phi_ks);
plot(subplot(4,4,7),t_vect.*0.001,I_ks)
title('I ks')
xlabel('s') 
ylabel('pA/pF')

I_pk =  Cmax_pk*(1+exp((25-V)./5.98)).^(-1).*(V-phi_k);
plot(subplot(4,4,8),t_vect.*0.001,I_pk)
title('I pk')
xlabel('s') 
ylabel('pA/pF')

I_to = Cmax_t0epi*G_r.*G_s.*(V-phi_k);
plot(subplot(4,4,9),t_vect.*0.001,I_to)
title('I to')
xlabel('s') 
ylabel('pA/pF')

I_cal = Cmax_cal.*G_d.*G_f.*G_fca*4*F^2.*V*(R*T)^(-1).*(C_ca.*exp(2*V*F*(R*T)^(-1))-0.341*c_ca0).*(exp(2*V*F*(R*T)^(-1))-1).^(-1);
plot(subplot(4,4,10),t_vect.*0.001,I_cal)
title('I cal')
xlabel('s') 
ylabel('pA/pF')

I_bca = Cmax_bca.*(V-phi_ca);
plot(subplot(4,4,11),t_vect.*0.001,I_bca)
title('I bca')
xlabel('s') 
ylabel('pA/pF')

I_pca = Cmax_pca.*C_ca.*(c_pca+C_ca).^(-1);
plot(subplot(4,4,12),t_vect.*0.001,I_pca)
title('I pca')
xlabel('s') 
ylabel('pA/pF')

I_leak = Imax_leak*(C_srca-C_ca);
plot(subplot(4,4,13),t_vect.*0.001,I_leak)
title('I leak')
xlabel('s') 
ylabel('mM/ms')

I_up = Imax_up*(1+c_up^2./C_ca.^2).^(-1);
plot(subplot(4,4,14),t_vect.*0.001,I_up)
title('I up')
xlabel('s') 
ylabel('mM/ms')

I_rel=Imax_rel.*G_d.*G_g.*(1+y_rel.*C_srca.^2.*(c_rel^2+C_srca.^2).^(-1));
plot(subplot(4,4,15),t_vect.*0.001,I_rel)
title('I rel')
xlabel('s') 
ylabel('mM/ms')

