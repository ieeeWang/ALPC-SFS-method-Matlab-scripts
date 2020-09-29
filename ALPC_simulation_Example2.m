clc; clear; close all;
% estimating latency of two systems with different time delay 
% for ALPC paper scenario 2

% declare the helper function
utils_1 = Utils1; 
utils_2 = Utils2; 

%% generate signals
% general settings
fs = 2000;
phase0 =  0.5*pi; %
ERPt = 0; % for simulation

tau1 =  0.051; % [sec]
tau2 =  0.021; % [sec]

w1 = [37, 43]; %Hz
w2 = [38, 46]; %Hz

[x1, y1_order2, t] = generateSignal_12sec(w1, 2, tau1, phase0, [], fs, ERPt);
[~, y1_order3, ~] = generateSignal_12sec(w1, 3, tau1, phase0, [], fs, ERPt);
y1 = y1_order2 + y1_order3;

[x2, y2_order2, ~] = generateSignal_12sec(w2, 2, tau2, phase0, [], fs, ERPt);
[~, y2_order3, ~] = generateSignal_12sec(w2, 3, tau2, phase0, [], fs, ERPt);
y2 = y2_order2 + y2_order3;

xx = x1 + x2;
yy = y1 + y2;

% add noise on yy
SNR_Y =  5; % dB
rng(1)
Y = awgn1(yy, SNR_Y,'measured','db');

%% check generated signals
% plot input & output signals
figure
 subplot(221);% Plot the EEG signal in the time domain.
    plot(t, xx)
    xlabel('t (s)');ylabel('Amp')  
    xlim([0 1])
    grid on
    title('X= x_1+x_2')
 subplot(222);
    [Y_f1, f] = utils_1.getfft_c(xx, fs); 
    plot(f, Y_f1)
    xlim([30 50]); ylim([0 1.2]) 
    title('X'); xlabel('f (Hz)');ylabel('|f(t)|')  
 subplot(223);% Plot the EEG signal in the time domain.
    plot(t, yy)
    xlabel('t (s)');ylabel('Amp')  
    xlim([0 1])
    grid on
    title('y_1+y_2 (\tau_1=51 ms, \tau_2=21 ms)')
 subplot(224);
    [Y_f1, f] = utils_1.getfft_c(yy, fs); 
    plot(f, Y_f1)
    xlim([0 140]); ylim([0 3])
    xlabel('f (Hz)');ylabel('|f(t)|') 
    title(['Y: ',num2str(SNR_Y), ' dB'])
    
%% output sets: 
% refer to 'demo_compute_high_order_outputs_compare2groups.m'
f1_2order = [6    74    80    86];
f1_3order = [31    37    43    49   111   117   123   129];
f2_2order = [8    76    84    92];
f2_3order = [30    38    46    54   114   122   130   138];

Yf = [f1_2order, f1_3order, f2_2order, f2_3order];
Yf_2order = [f1_2order, f2_2order];
Yf_3order = [f1_3order, f2_3order];
Yf_1 = [f1_2order, f1_3order];
Yf_2 = [f2_2order, f2_3order];


%% plot one original subsystem
% % f_output = Yf_1;
% f_output = Yf_2;
% % f_output = Yf_2order;
% % f_output = Yf_3order;
% 
% [f_output, ~] = sort(f_output,'ascend');% ranking 'f_output'
% plt = 1;
% [min_MSE, tau_est, ~] = apprent_latency_cycles(f_output, Y, [], fs, plt, -1*ERPt);

%% plot two original subsystems in one figure
plt = 0;
[min_MSE1, tau_est1, pltOut1] = utils_1.apprent_latency_cycles(Yf_1, Y, [], fs, plt, -1*ERPt);
[min_MSE2, tau_est2, pltOut2] = utils_1.apprent_latency_cycles(Yf_2, Y, [], fs, plt, -1*ERPt);

utils_2.plot2fitlines_inOne(pltOut1, pltOut2)

%% plot MSPC in one figure: given same input and output, MSPC = ALPC
plt = 0;
[min_MSE1, tau_est1, pltOut1_mspc] = utils_1.apprent_latency_cycles(Yf_2order, Y, [], fs, plt, -1*ERPt);
[min_MSE2, tau_est2, pltOut2_mspc] = utils_1.apprent_latency_cycles(Yf_3order, Y, [], fs, plt, -1*ERPt);

utils_2.plot2fitlines_inOne(pltOut1_mspc, pltOut2_mspc)

%% ALPC + SFS
% (1) get phase angle on each output f
f_output = Yf;
[f_output, ~] = sort(f_output,'ascend');% ranking 'f_output'
P_compx = []; % get FFT complex value at 'f_input'
for i=1:length(f_output)
    P_compx(i) = utils_1.getfft_compx(Y, fs, f_output(i)); 
end
angle_f_out = angle(P_compx); % angles lie between [-pi, +pi]
f_angle = [f_output(:), angle_f_out(:)];

% (2) ALPC + forward subset selection (SFS)
% Option A: SFS starting from 1 single f
f1 = 6

TC = -1*ERPt; % see my paper for TC
tmax = 0.15;
[output2] =  utils_1.getALPC_SFS1_TC(f1, f_angle, TC, tmax,1);
f_sel = output2.fsel;            
        
% use the selected subset
n0 = 11;% the step where MSE starts to ascend
f_subset = f_sel(1: n0+1);
[min_MSE, tau_est, pltOut1_alpc]  = utils_1.apprent_latency_cycles(f_subset, yy, [], fs, 1, -1*ERPt);
fprintf('sub-system 1...')
min_MSE
tau_est

% (1) use the remaining as subset2
f_subset = f_sel(n0+2: end);
[min_MSE, tau_est, pltOut2_alpc]  = utils_1.apprent_latency_cycles(f_subset, yy, [], fs, 1, -1*ERPt);
fprintf('sub-system 2...')
min_MSE
tau_est


%% plot ALPC+SFS in one figure
utils_2.plot2fitlines_inOne(pltOut1_alpc, pltOut2_alpc)

%% ALPC paper figure plot
utils_2.plot_MSPC_ALPC_inOne(pltOut1_alpc, pltOut2_alpc, pltOut1_mspc, pltOut2_mspc)














