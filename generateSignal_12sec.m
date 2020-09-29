function [xx, yy, t] = generateSignal_12sec (w, order, tau, phase0, SNR_Y, fs, ERPt)
%  generate sianal of given frequencies and orders with 12 sec.
% INPUT
%     w - given frequencies
%     order - system order
%     tau - system delay
%     phase0 - initial pahse
%     SNR_Y - SNR
%     fs - sampling ratio
%     ERPt - stimuli delay, default = 0 
% OUTPUT
%     xx - sum of all inputs
%     yy - it MUST be a 12-s time series, the function 'getfft_compx.m' is
%     based on 1/12 f resolustion. Otherwise change this function.
%     t - time labels for plotting

%% [input  signals] 
tt = 12; % s
deltaT = 1* ERPt*fs; % make signal go ahead << correct for real stimuli (0.3s)

T = 1/fs;             % Sampling period       
L = fs* tt;             % Length of signal
t = (0:L-1)*T;        % Time vector
t_real = ((0:L-1)+ deltaT)*T;  % time of input signal


% input signals
xx=0;
for i =1:length(w) 
    xx = xx + sin(2*pi*w(i)*t_real + phase0); % with deltaT delay
end

% get FFT complex value at 'w'
% decare the helper function
Utils = Utils1; 

P_compx_w = [];
for i=1:length(w)
    P_compx_w(i) = Utils.getfft_compx(xx, fs, w(i)); %N = 33 %Hz
end
% check if the 'angle_f_in' == 0?
angle_f_in = angle(P_compx_w);


%% output signals 
t_tmp = t_real - tau; 
yi=0;
for i =1:length(w) 
    yi = yi + sin(2*pi*w(i)* t_tmp + phase0); % with deltaT delay
end
% yy = yi;
yy = yi.^order;


if ~isempty(SNR_Y)
    yy = awgn1(yy, SNR_Y,'measured','db');
end

% check sine on the minus time, which is still continuous
% figure
%     plot(t_tmp, sin(2*pi*7*t_tmp))
%     grid on
%     ylim([-1.5 1.5])

