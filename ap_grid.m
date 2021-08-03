%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Simulate the Aliev-Paniflov Model for a grid of cells          %
%                              Forward Euler                              %
%                                                                         %
%                          Program Name: ap_grid                          %
%                                                                         %
%                         History: Created 06/07                          %
%                                  Update 07/07                           %
%                                                                         %
%                                                                         %
%                       Author: Alix Vanpoperinghe                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all
%creates folders in location of matlab file
%mkdir frames_ap
%mkdir movie_ap

%% PARAMETERS
len_grid = 100;
lx = 100;
dx = lx/len_grid;

e0 = 0.002;
k = 8;
b = 0.05;                                       % oscillation threshold (positive = stable)
mu_1 = 0.2;
mu_2 = 0.3;
epsilon = @(u,v) e0 + (mu_1*v)./(u + mu_2);     % weighting factor

diff = 0.1;                                     % diffusion coefficient

iapp_stim = 0.1;                                % stimulus strength
iapp_width = 2;                                 % width of cells being stimulated

t_burst_start = 0;                              % start time of Iapp burst
t_burst_fin = 2;                                % end time of Iapp burst

%% POSITION
% select cell to investigate
cell_x = 60;
cell_y = 50;
pos = [cell_x cell_y];

%% AP
du = @(u,v,iapp,lap) - k*u.*(u-b).*(u-1) -u.*v + diff*lap + iapp;       % fast
dv = @(u,v) epsilon(u,v).*(-v - k*u.*(u-b-1));                               % slow

%% INTEGRATION SETUP
t_start = 0;                    % start time
dt = 0.01;                      % integration step size
t_step = 0.5;                   % recording step size

t_fin = 100;                    % end time
N = t_fin/dt;                   % number of grid points
M = t_fin/t_step;               % number of frames

t = linspace(t_start, t_fin, N);% storing time values


%% STEADY STATE
t_reach = 1000;
u = zeros(len_grid,len_grid);
v = zeros(len_grid,len_grid);
SS = steadyState(u,du,v,dv,dt,t_reach/dt);
tspan = dt:dt:t_reach;



%% INITIAL CONDITIONS
u = SS(1,end)*ones(len_grid,len_grid);
v = SS(2,end)*ones(len_grid,len_grid);
ep = zeros(len_grid,len_grid);
x = 1:dx:lx;
y = 1:dx:lx;
iapp = zeros(len_grid,len_grid);
S = zeros(len_grid,len_grid,2);
cell = zeros(1,N);


%% OUTPUT
% change to 1 to select action(s)
draw_fig = 1;                       %do you want figures to display in matlab?
save_im = 0;                        %frames
save_mov = 0;                       %movie
plot_ts = 0;                        %cell time series
import_csv = 0;                     %import iapp from csv file

% name movie here

names = randi(100);
movieFileName = sprintf('movie_%03d', names);
%movieFileName = 'rand_99_long';
fullFileName2 = fullfile(pwd,'movie_ap',movieFileName);
vid = VideoWriter(fullFileName2, 'MPEG-4');

%% DISK
% creates a disk of 1s in an array of zeros
mid = len_grid/2;

p = mid;
q = mid;
radius = 20;

[X, Y] = meshgrid(-p+1:len_grid-p, -q+1:len_grid-q);
iapp_disk = X.^2 + Y.^2 <= radius^2;

%% IAPP FROM CSV
if (import_csv == 1)
    %input file name here. Must be in same location as this matlab file.
    T = readtable('myfile.csv');
    [j,k] = size(T);
    if (j==k && j==len_grid)            % assuming here input is only 0s and 1s
        iapp_csv(:,:) = T{:,:};
    else
        disp('Incorrect dimensions in csv file');
    end
end

%% IAPP RANDOM
B = zeros(len_grid,len_grid);
ran = rand(len_grid);
for len = 1:len_grid
    for wid = 1:len_grid
        if (ran(len,wid) >= 0.99)
            B(len,wid) = 1;
        end
    end
end
iapp_rand = B;

%% FORWARD EULER
open(vid);
for it = 1:N-1
        if(it*dt>=t_burst_start && it*dt<t_burst_fin)
            iapp(:,1:3) = iapp_stim;
            %iapp = iapp_stim*iapp_disk;
            %iapp = iapp_stim*iapp_csv;
            %iapp = iapp_stim*iapp_rand;
        elseif (it*dt>=110 && it*dt<=112)
            iapp(1:3,:) = iapp_stim;
        else
            iapp = zeros(len_grid,len_grid);
        end
        S = eulerStep(u, du, v, dv, iapp, dx, dt);
        u = S(:,:,1);
        v = S(:,:,2);
        cell(it) = S(pos(1), pos(2), 1);
        
        
        
        
        if (rem(it*dt, t_step) == 0 && draw_fig == 1)
            h = surf(x, y, S(:,:,1), 'EdgeColor', 'none');
            
            %formatting figure
            view(2);
            shading interp
            grid off
            caxis([0 1]);
            colorbar
            set(gca, 'XTick', [], 'YTick', []);
            txt = ['t = ', num2str(it*dt,'%3.1f')];
            title(txt);
            drawnow;
            
            %save frames
            if (save_im == 1)
                thisFrame = getframe(gcf);
                NumFrame = it*dt/t_step;
                baseFileName = sprintf('frame_%03d.png',NumFrame);
                fullFileName = fullfile(pwd,'frames_ap',baseFileName);
                saveas(gcf,fullFileName);
            end
            
            %convert frames into video            
            if (save_mov == 1)
                writeVideo(vid,thisFrame);
            end
        end
end

close(vid);

%% PLOTTING
if (plot_ts == 1)
    plot(t,cell(1,:));
    xlabel('Time');ylabel('Voltage');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

