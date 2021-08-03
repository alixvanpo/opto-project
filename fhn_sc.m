%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Simulate the FitHugh-Nagumo Model for a single cell           %
%                              Forward Euler                              %
%                                                                         %
%                          Program Name: fhn_sc                           %
%                                                                         %
%                       History: Created 18/06/2021                       %
%                                                                         %
%                       Author: Alix Vanpoperinghe                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PARAMETERS
epsilon = 0.01;                  % "abruptness" of excitation = 0.01
a = 0.15;                       % threshold of excitation = 0.5
b = 0.5;                        % 0.5 (0.32 for longer plateau)
c = 1;                          % 1
d = 0;                          % 0

g = 0.01;       
p = 0.7;        
q = 0.7;        

t_burst_start = 0;              % start time of Iapp burst
t_burst_fin = 2;                % end time of Iapp burst

dv = @(v,w,iapp) -v*(v-a)*(v-1) - w + iapp;               % fast
dw = @(v,w) epsilon*(b*v + d - c*w);                      % slow

%% INTEGRATION SETUP
t_start = 0;                    % start time
dt = 0.01;                      % integration step size
t_step = 0.5;                   % recording step size
t_fin = 100;                    % end time
N = t_fin/dt;                   % number of grid points

t = linspace(t_start, t_fin, N);% storing time values

%% STEADY STATE
t_reach = 1000;
v = zeros(N,1);
w = zeros(N,1);
SS = steadyState_sc(v,dv,w,dw,dt,t_reach/dt);
tspan = dt:dt:t_reach;


%% INITIAL CONDITIONS
v = SS(1,end)*ones(N,1);
w = SS(2,end)*ones(N,1);
iapp = zeros(N,1);

%% FORWARD EULER
for it = 1:N-1
        if(it*dt>=t_burst_start && it*dt<t_burst_fin)
            iapp(it) = 0.5;
        else 
            iapp(it) = 0;
        end
        v(it+1) = v(it) + dt*dv(v(it),w(it),iapp(it));
        w(it+1) = w(it) + dt*dw(v(it),w(it));
	end


%% PLOTTING
plot(t, v, 'r')
xlabel('Time'),ylabel('Voltage')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

