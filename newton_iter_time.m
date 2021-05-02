function [V_nou] = newton_iter_time(V_ant,gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,gnou_r,gnou_s,gnou_d,gnou_f,gnou_k1,gnou_fca,...
                    cnou_na,cnou_k,cnou_ca,dt,n)

constants

phi_k =  phi_ion(R,T,F,z_k,c_k0,cnou_k);
phi_na = phi_ion(R,T,F,z_na,c_na0,cnou_na);
phi_ca = phi_ion(R,T,F,z_ca,c_ca0,cnou_ca);
phi_ks = R*T/F*log((c_k0+p_kna*c_na0)*(cnou_k+p_kna*cnou_na)^(-1));

% Initialize currents with potential as x variable
I_na =  @(x) Cmax_na*gnou_m^3*gnou_h*gnou_j*(x-phi_na);
I_bna = @(x) Cmax_bna*(x-phi_na);
I_nak = @(x) Imax_nak*(c_k0*cnou_na)*((cnou_na+c_nak)*(c_k0+c_kna)*(1+0.1245*exp(-0.1*x*F/(R*T))+0.0353*exp(-x*F/(R*T))))^(-1);
I_naca = @(x) Imax_naca*(exp(y*x*F/(R*T))*cnou_na^3*c_ca0-exp((y-1)*x*F/(R*T))*c_na0^3*cnou_ca*y_naca)...
            *((c_naca^3+c_na0^3)*(c_cana+c_ca0)*(1+k_naca*exp((y-1)*x*F/(R*T))))^(-1);
I_k1 = @(x) Cmax_k1*gnou_k1*(c_k0/5.4)^(1/2)*(x-phi_k);
I_kr = @(x) Cmax_kr*gnou_xr1*gnou_xr2*(c_k0/5.4)^(1/2)*(x-phi_k);
I_ks = @(x) Cmax_ksepi*gnou_xs^2*(x-phi_ks);
I_pk = @(x) Cmax_pk*(1+exp((25-x)/5.98))^(-1)*(x-phi_k);
I_t0 = @(x) Cmax_t0epi*gnou_r*gnou_s*(x-phi_k);
I_cal = @(x) Cmax_cal*gnou_d*gnou_f*gnou_fca*4*F^2*x*(R*T)^(-1)*(cnou_ca*exp(2*x*F*(R*T)^(-1))-0.341*c_ca0)*(exp(2*x*F*(R*T)^(-1))-1)^(-1);
I_bca = @(x) Cmax_bca*(x-phi_ca);
I_pca = Cmax_pca*cnou_ca*(c_pca+cnou_ca)^(-1);

% Initialize potential function
f = @(x)(x-V_ant+(I_na(x)+I_bna(x)+I_nak(x)+I_naca(x)+I_k1(x)+I_kr(x)+I_ks(x)+I_pk(x)+I_t0(x)+I_cal(x)+I_bca(x)+I_pca- stim(n,dt))*dt);

% Derivatives of current functions
dx_I_na = Cmax_na*gnou_m^3*gnou_h*gnou_j;
dx_I_bna = Cmax_bna;
dx_I_nak = @(x)Imax_nak*c_k0*cnou_na*((cnou_na+c_nak)*(c_k0+c_kna))^(-1)...
            *(-(F*(R*T)^(-1)*(-0.0353*exp(-F*(R*T)^(-1)*x)-0.01245*exp(-0.1*F*(R*T)^(-1)*x)))/(1+0.0353*exp(-F*(R*T)^(-1)*x)+0.1245*exp(-0.1*F*(R*T)^(-1)*x))^2);
dx_I_naca = @(x)(Imax_naca*F*(R*T)^(-1)*exp((y-1)*x*F*(R*T)^(-1)))...
                 *(cnou_na^3*c_ca0*(k_naca*exp(y*x*F*(R*T)^(-1))+y*exp(x*F*(R*T)^(-1)))-c_na0^3*cnou_ca*y_naca*(y-1))...
                 *((c_naca^3+c_na0^3)*(c_cana+c_ca0))^(-1)*(1+k_naca*exp((y-1)*x*F*(R*T)^(-1)))^(-2);
dx_I_k1 = Cmax_k1*gnou_k1*(c_k0/5.4)^(1/2);
dx_I_kr = Cmax_kr*gnou_xr1*gnou_xr2*(c_k0/5.4)^(1/2);
dx_I_ks = Cmax_ksepi*gnou_xs^2;
dx_I_pk = @(x)(Cmax_pk*exp(0.167224*x)*(65.4052-10.9373*phi_k+exp(0.167224*x)+10.9373*x))*(65.4052+exp(0.167224*x))^(-2);
dx_I_t0 =  Cmax_t0epi*gnou_r*gnou_s;
dx_I_cal = @(x)Cmax_cal*gnou_d*gnou_f*gnou_fca*4*F^2*(R*T)^(-1)*(cnou_ca*exp(2*F*(R*T)^(-1)*x)*(-2*F*(R*T)^(-1)*x+exp(2*F*(R*T)^(-1)*x)-1)+c_ca0*(exp(2*F*(R*T)^(-1)*x)*(0.341*2*F*(R*T)^(-1)*x-0.341)+0.341))...
                *(exp(2*F*(R*T)^(-1)*x)-1)^(-2);
dx_I_bca = Cmax_bca;
dx_I_pca = 0;

% Initialize derivative of potential function
df = @(x)(1+dt*(dx_I_na + dx_I_bna + dx_I_nak(x) + dx_I_naca(x) + dx_I_k1 + dx_I_kr + dx_I_ks + dx_I_pk(x) + dx_I_t0 + dx_I_cal(x) + dx_I_bca + dx_I_pca));

% Initial values
x0 = V_ant; f0 = f(x0);
df0 = df(x0); 

% Iteration
x1 = x0 - f0/df0; 
V_nou = x1;
% V_nou=V_ant-f0;


