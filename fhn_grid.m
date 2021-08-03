%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Simulate the FitHugh-Nagumo Model for a grid of cells          %
%                              Forward Euler                              %
%                                                                         %
%                         Program Name: fhn_grid                          %
%                                                                         %
%                          History: Created 22/06                         %
%                                   Disk 06/07                            %
%                                   Time Series 07/07                     %
%                                   Iapp shapes 07/07                     %
%                                   Steady State 09/07                    %
%                                                                         %
%                                                                         %
%                       Author: Alix Vanpoperinghe                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
%creates folders in location of matlab file
mkdir frames_fhn
mkdir movie_fhn

%% PARAMETERS
epsilon = 0.01;                  % "abruptness" of excitation
a = 0.15;                        % threshold of excitation
b = 0.5;
c = 1;
d = 0.5;                           % making this nonzero removes oscillations!! what are these oscs?
diff = 2;                        % diffusion coefficient

len_grid = 100;                   % number of points in length
lx = 100;                         % dimensional length
dx = lx/len_grid;

t_burst_start = 0;               % start time of Iapp burst
t_burst_fin = 3;                 % end time of Iapp burst

iapp_width = 1:20;                % stimulus width of cells
iapp_height = 1:20;
iapp_stim = 3;                   % stimulus strength


%% POSITION
% select cell to investigate
cell_x = 16;
cell_y = 16;
pos = [cell_x cell_y];

%% FHN
%FitzHugh-Nagumo formalism
dv = @(v,w,iapp,lap) -v.*(v-a).*(v-1) - w + diff*lap + iapp;              % fast
dw = @(v,w) epsilon*(b*v + d - c*w);                                      % slow

%% INTEGRATION SETUP
t_start = 0;                    % start time
dt = 0.01;                      % integration step size
t_step = 0.1;                   % recording step size


t_fin = 50;                     % end time
N = t_fin/dt;                   % number of grid points
M = t_fin/t_step;               % number of frames

t = linspace(t_start, t_fin, N);% storing time values

%% STEADY STATE
t_reach = 1000;
v = zeros(len_grid,len_grid);
w = zeros(len_grid,len_grid);
SS = steadyState(v,dv,w,dw,dt,t_reach/dt);
tspan = dt:dt:t_reach;


%% INITIAL CONDITIONS
v = SS(1,end)*ones(len_grid,len_grid);
w = SS(2,end)*ones(len_grid,len_grid);
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

%names = randi(100);
%movieFileName = sprintf('disk_%03d', names);
movieFileName = 'disk_works';
fullFileName2 = fullfile(pwd,'movie_fhn',movieFileName);
vid = VideoWriter(fullFileName2, 'MPEG-4');


%% DISK
% creates a disk of 1s in an array of zeros
mid = len_grid/2;

p = mid;                    %(p,q) is the location of the center of the disk
q = mid;
radius = 20;

[X, Y] = meshgrid(-p+1:len_grid-p, -q+1:len_grid-q);
iapp_disk = X.^2 + Y.^2 <= radius^2;

%% IAPP FROM CSV
if (import_csv == 1)
    T = readtable('myfile.csv');
    [j,k] = size(T);
    if (j==k && j==len_grid)            % assuming here input is only 0s and 1s
        iapp_csv(:,:) = T{:,:};
    end
end

%% IAPP RANDOM
B = zeros(len_grid,len_grid);
constraint = 0.85;                  %percentage of grid that will be stimulated
ran = rand(len_grid);
for len = 1:len_grid
    for wid = 1:len_grid
        if (ran(len,wid) >= constraint)
            B(len,wid) = 1;
        end
    end
end
iapp_rand = B;

%% FORWARD EULER
open(vid);
for it = 1:N-1
        if(it*dt>=t_burst_start && it*dt<t_burst_fin)
            %iapp(:,iapp_width) = iapp_stim;
            iapp = iapp_stim*iapp_disk;
            %iapp = iapp_stim*iapp_csv;
            %iapp = iapp_stim*iapp_rand;
        %elseif (it*dt>=t_burst_fin && it*dt<(t_burst_fin+2)) 
            %iapp = zeros(len_grid,len_grid);
            %iapp = -iapp_stim*iapp_disk;
            %iapp = -SS(1,end)*ones(len_grid,len_grid);
            %iapp = 0.4*ones(len_grid,len_grid) - 0.4*iapp_disk;
            %iapp = 0.2*ones(len_grid,len_grid);
        else
            iapp = 0.5*ones(len_grid,len_grid);
            %iapp = zeros(len_grid,len_grid);
        end
        
        %sets current state of grid
        S(:,:,1) = v;
        S(:,:,2) = w;
        
        %calculates next state (whole grid) and updates v,w
        S_next = eulerStep(v, dv, w, dw, iapp, dx, dt);
        v = S_next(:,:,1);
        w = S_next(:,:,2);
        
        %outputs current state of cell
        cell(it) = S(pos(1), pos(2), 1);
        
        
        
        if (rem(it*dt, t_step) == 0 && draw_fig == 1)
            h = surf(x, y, S(:,:,1), 'EdgeColor', 'none');
            
            %formatting figure
            view(2);
            shading interp
            grid off
            caxis([-0.4 1]);
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
                fullFileName = fullfile(pwd,'frames_fhn',baseFileName);
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

%plot(tspan,SS(1,:))
