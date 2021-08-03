%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Simulate the Aliev-Panfilov Model for a single cell           %
%                              Forward Euler                              %
%                                                                         %
%                         Program Name: aliev_sc                          %
%                                                                         %
%                       History: Created 22/06/2021                       %
%                                                                         %
%                       Author: Alix Vanpoperinghe                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PARAMETERS
t_burst_start = 0;              % start time of Iapp burst
t_burst_fin = 2;                % end time of Iapp burst
%iapp = 0.03;                   % applied stimulus

e0 = 0.002;
k = 8;
b = 0.05;                       % oscillation threshold (positive = stable)
mu_1 = 0.2;
mu_2 = 0.3;
epsilon = @(u,v) e0 + (mu_1*v)/(u + mu_2);            % weighting factor

du = @(u,v,iapp) - k*u*(u-b)*(u-1) -u*v + iapp;       % fast
dv = @(u,v) epsilon(u,v)*(-v - k*u*(u-b-1));          % slow

%% INTEGRATION SETUP
t_start = 0;                    % start time
dt = 0.01;                      % integration step size
t_step = 0.5;                   % recording step size
t_fin = 100;                    % end time
N = t_fin/dt;                   % number of grid points

tau = linspace(t_start, t_fin, N);% storing time values

%% STEADY STATE
t_reach = 1000;
u = zeros(N,1);
v = zeros(N,1);
SS = steadyState_sc(u,du,v,dv,dt,t_reach/dt);
tspan = dt:dt:t_reach;


%% INITIAL CONDITIONS
u = SS(1,end)*ones(N,1);
v = SS(2,end)*ones(N,1);
iapp = zeros(N,1);

%% FORWARD EULER
for it = 1:N-1
    if (it*dt>=t_burst_start && it*dt<t_burst_fin)
            iapp(it) = 0.5;
        else 
            iapp(it) = 0;
    end
        u(it+1) = u(it) + dt*du(u(it),v(it),iapp(it));
        v(it+1) = v(it) + dt*dv(u(it),v(it));
	end


%% PLOTTING
%vol = 100*u - 80;               % 80 - 87.3 mV (scaling)
%t = 7*tau;                      % (scaling)
plot(tau, u, 'r')
xlabel('Time (ms)'),ylabel('Voltage (mV)')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

