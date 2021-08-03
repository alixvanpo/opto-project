

%  plotting comparative and investigative time series of opto file   %



%% COMPARING STIMULI
close all
clear
clc

t = 0.01:0.01:100;
tau = 4*t;


%LIGHT STIMULATION
v = fhn_comp(0,0.7);                        % 0.5 mW/mm2, 470 nm for 10 ms (faster rise = 1)
% v = ap_comp(0,0.7);
vol = 100*v - 80;
h(1) = plot(tau, vol, 'DisplayName', 'OPTICAL');


hold on

% ELECTRICAL STIMULATION
v = fhn_comp(0.5,0);                        % 8pA/pF for 5ms
% v = ap_comp(0.5,0);
vol = 100*v - 80;
h(2) = plot(tau, vol, 'DisplayName', 'ELECTRICAL');

title('Comparing Stimuli');
legend(h);
xlabel('Time'),ylabel('Voltage')

hold off;



%% TRACKING STATE OF CHANNEL
close all
t = 0.01:0.01:100;
tau = 4*t;
[v,w,C1,O1,O2,C2] = ap_comp(0,0.4);        % 0.5 mW/mm2, 470 nm for 10 ms (faster rise = 1)
tr(1) = plot(tau, C1(:), 'DisplayName', 'closed sensitized');

hold on
tr(2) = plot(tau, O1(:), 'DisplayName', 'excited');
tr(3) = plot(tau, O2(:), 'DisplayName', 'open');
tr(4) = plot(tau, C2(:), 'DisplayName', 'closed desensitized');


title('State of ChR2');
legend(tr);
xlabel('Time'),ylabel('State Fraction')

hold off;












