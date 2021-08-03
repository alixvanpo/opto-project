%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Aliev-Panfilov Single Cell Model with Optogenetics            %
%                              Forward Euler                              %
%                                                                         %
%                        Program Name: ap_sc_opto                         %
%                                                                         %
%                         History: Created 26/07                          %
%                                                                         %
%                                                                         %
%                       Author: Alix Vanpoperinghe                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
% close all
% clc

%% AP FULL OPTOGENETIC MODEL
% AP PARAMETERS
e0 = 0.002;
k = 8;
b = 0.05;                       % oscillation threshold (positive = stable)
mu_1 = 0.2;
mu_2 = 0.3;
epsilon = @(u,v) e0 + (mu_1*v)/(u + mu_2);            % weighting factor

% ChR2 CURRENT
g_chr = 0.75;                                        % max conductance of ChR (val 0.4-0.75 mS/cm2)
exp_lvl = 1;                                         % expression level of ChR2, default is 1
gamma = 0.1;                                         % ratio O1/O2 (val 0.05-0.1)
G = @(u) (10.6408 - 14.6408*exp(-u/42.7671))/u;      % voltage rectification function (NOT SURE ABOUT DIMENSIONS HERE)
%vol = @(u) 100*u - 80;
E_chr = 0;                                           % reversal potential of ChR (val -5 - +10 mV)

ichr = @(u,O1,O2) exp_lvl*g_chr*(O1 + gamma*O2)*G(u)*(u-E_chr);     % ChR current

% FOUR-STATE MODEL PARAMETERS AND DEs
%k1 = @(F,p) ep1*F*p;
k1 = 0.369;                          % light-sensitive rate constant C1 to O1
%k2 = @(F,p) ep2*F*p;
k2 = 0.369;                          % light-sensitive rate constant C2 to O2
%F = sig_ret * blabla
%I = 10;                             % irradiance (values 0-10)
%lambda = 470;                       % (nm) wavelength of max absorption
%w_loss = 0.77;                      % scaling factor for loss of photons
%F = 0.00006*I*(lambda/w_loss);

% rate constants
Gd1 = 0.084;                        % O1 to C1   (val 0.084-0.11)
Gd2 = 0.1254;                       % O2 to C2   (val 0.025-0.1254)
Gr = 0.004;                         % C2 to C1   (val 0.0004-0.0040)
e12 = 0.011;                        % O1 to O2   (val 0-0.03)
e21 = 0.008;                        % O2 to O1   (val 0.008-0.015)



dC1 = @(O1,C1,C2) Gr*C2 + Gd1*O1 - k1*C1;
dO1 = @(O1,O2,C1) k1*C1 - (Gd1 + e12)*O1 + e21*O2;
dO2 = @(O1,O2,C2) k2*C2 - (Gd2 + e21)*O2 + e12*O1;
dC2 = @(O2,C2) Gd2*O2 - (k2 + Gr)*C2;


% FHN DEs
Cm = 60;                           % (pF) 56.6 +/- 9.9 - Hotka et al. (2014)
Rm = 1/Cm;

du = @(u,v,iapp,O1,O2) - k*u*(u-b)*(u-1) -u*v + iapp - Rm*ichr(u,O1,O2);       % fast
dv = @(u,v) epsilon(u,v)*(-v - k*u*(u-b-1));                                   % slow

%% INTEGRATION SETUP
t_start = 0;                    % start time
dt = 0.01;                      % integration step size



t_fin = 100;                    % end time
N = t_fin/dt;                   % number of grid points

t = linspace(t_start, t_fin, N);% storing time values


t_burst_start = 0;              % start time of Iapp burst
t_burst_fin = 2;                % end time of Iapp burst



%% STEADY STATE?
 t_reach = 1000;
 delta = 0.01;
 M = t_reach/delta;

u = ones(M,1);
v = ones(M,1);
for s = 1:M-1
    u(s+1) = u(s) + delta*du(u(s),v(s),0,0,0);
    v(s+1) = v(s) + delta*dv(u(s),v(s));
end
tspan = delta:delta:t_reach;

%% INITIAL CONDITIONS
u = u(end)*ones(N,1);
v = v(end)*ones(N,1);

C1 = ones(N,1);                     % assumed 100% of ChRs are in the C1 state
O1 = zeros(N,1);
O2 = zeros(N,1);
C2 = zeros(N,1);
iapp = zeros(N,1);

I_CHR2 = zeros(N,1);                % current tracker

%% FORWARD EULER
for it = 1:N-1
        u(it+1) = u(it) + dt*du(u(it),v(it),iapp(it),O1(it),O2(it));
        v(it+1) = v(it) + dt*dv(u(it),v(it));
        C1(it+1) = C1(it) + dt*dC1(O1(it),C1(it),C2(it));
        O1(it+1) = O1(it) + dt*dO1(O1(it),O2(it),C1(it));
        O2(it+1) = O2(it) + dt*dO2(O1(it),O2(it),C2(it));
        C2(it+1) = C2(it) + dt*dC2(O2(it),C2(it));
        
        
        I_CHR2(it+1) = ichr(u(it+1),O1(it+1),O2(it+1));
	end


%% PLOTTING
vol = 100*u - 80;               % 80 - 87.3 mV (scaling)
tau = 4*t;                      % (scaling)
plot(tau, vol, 'r')
xlabel('Time'),ylabel('Voltage')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

