clc; clear; close all;
% ALPC-SFS paper Example 1: 
% Two pure second-order subsystems with non-overlapping inputs.

% declare the helper function
Utils = Utils1; 

%% generate X and Y
order = 2;
fs = 2000;
phase0 =  0.5*pi; %
w1 = [17, 21, 27]; %Hz
tau1 = 0.051; % [sec]

% ERPt = 0.3; % EEG trials start at 0.3s (stimulus delay = 0.3s) in real experiment (Fig. 5)
ERPt = 0; % for simulation Example 1 without stimulus delay (Fig. 4)
[xx1, yy1, t] = generateSignal_12sec(w1, order, tau1, phase0, [], fs, ERPt);

w2 = [41, 49]; %Hz
tau2 = 0.021; % [sec]
[xx2, yy2, ~] = generateSignal_12sec(w2, order, tau2, phase0, [], fs, ERPt);

xx = xx1 + xx2;
yy0 = yy1 + yy2;

% add AWGN noise
SNR_Y =  5; % dB
rng(1)
yy = awgn1(yy0, SNR_Y,'measured','db');
    
%% plot input & output signals
% (1) visulize input(sound stimuli)
figure
fmax = 100; % Hz
 subplot(221);% Plot the EEG signal in the time domain.
    plot(t, xx)
    xlabel('t (s)');ylabel('Amp')  
    xlim([0 1])
    grid on
    title('X= x_1+x_2')
 subplot(222);
    [Y_f1, f] = Utils.getfft_c(xx, fs); 
    plot(f, Y_f1)
    title('X'); xlabel('f (Hz)');ylabel('|f(t)|') 
    xlim([0 fmax]) 
 subplot(223);% Plot the EEG signal in the time domain.
    plot(t, yy0)
    xlabel('t (s)');ylabel('Amp')  
    xlim([0 1])
    grid on
    title('y_1+y_2 (\tau_1=51 ms, \tau_2=21 ms)')
 subplot(224);
    [Y_f1, f] = Utils.getfft_c(yy, fs); 
    plot(f, Y_f1)
    xlabel('f (Hz)');ylabel('|f(t)|') 
    title(['Y: ',num2str(SNR_Y), ' dB'])
    xlim([0 fmax]) 
        
%% (A) MSPC 
w = [w1, w2];

X = reshape(xx,[fs,12]);
Y = reshape(yy,[fs,12]);

MS = fft(X);
D2_delay = fft(Y);
 
% (1) seperately use [w1] and [w2] as input to get the outputs
[~,Cp1,~,~,f_sigma1] = Utils.MSPC_2b(MS,D2_delay,w1+1); 
[~,Cp2,~,~,f_sigma2] = Utils.MSPC_2b(MS,D2_delay,w2+1); 

C_p_2 = [Cp1; Cp2];
f_sigma_2 = [f_sigma1;f_sigma2]-1;

% (2) compute delay
mytitle = 'MSPC';
[min_MSE, tau_est] = Utils.MSPC_latency(f_sigma_2, C_p_2, 1, mytitle)

%% (B) ALPC without SFS
f_output = f_sigma_2;
% ranking 'f_output'
[f_output, ~] = sort(f_output,'ascend');
plt = 1;
[min_MSE, tau_est, ~] = Utils.apprent_latency_cycles(f_output, yy, [], fs, plt, -1*ERPt)

%% plot SNR spectrum
% sn = 12;
% figure   
%     [Y_f1, f] = getfft_c(yy, fs);   
%     SNR = getSNR_spectrum_c(f, Y_f1, 12, 0.5);
%     SNRdB = 10*log10(SNR);
%     subplot(311)
%         stem(SNR) % power ratio
%         grid on
%         xlabel('f (Hz)','Fontsize',sn);ylabel('SNR','Fontsize',sn) 
%         title(['SNR (target/neighbors) on integer f: ',num2str(SNR_Y), ' dB'])
%         xlim([0 fmax])
%     subplot(312)
%         stem(f_sigma1-1, ones(size(f_sigma1)),'r')
%         hold on
%          stem(f_sigma2-1, ones(size(f_sigma2)),'b')
%         legend('sys 1', 'sys 2')
%         grid on
%         xlabel('f (Hz)','Fontsize',sn);ylabel('output','Fontsize',sn) 
%         xlim([0 fmax])     
%     subplot(313)
%         targetdB = SNRdB(f_output);
% %         stem(f_output,targetdB)% dB
%         stem(SNRdB)% dB
%         grid on
%         xlabel('f (Hz)','Fontsize',sn);ylabel('SNR (dB)','Fontsize',sn) 
%         xlim([0 fmax])
%         ylim([0 30])
%  
% % SNR (dB) compaired with neighbering f bins        
% SNRdB_m = mean(targetdB)

%% get phase angle on each output f
P_compx = []; % get FFT complex value at 'f_input'
P_compx0 = []; 
for i=1:length(f_output)
    P_compx(i) = Utils.getfft_compx(yy, fs, f_output(i)); 
    P_compx0(i) = Utils.getfft_compx(yy0, fs, f_output(i)); 
end
angle_f_out = angle(P_compx); % angles lie between [-pi, +pi]
angle_f_out0 = angle(P_compx0); % angles lie between [-pi, +pi]
f_angle = [f_output(:), angle_f_out(:)];

% plot phase 
% figure
%     p0_vec = P_compx0./abs(P_compx0);
%     p_vec  = P_compx ./abs(P_compx);
%     PEL = abs(p_vec-p0_vec); % [0-2]
%     normError = pi*PEL/2;
%     PError_m = mean(normError); % mean PE on all target f
%     subplot(211)
%         stem(f_output, angle_f_out,'o')
%         hold on
%         stem(f_output, angle_f_out0,'*')
%         ylim([-pi pi])
%         grid on
%         legend('phase with noise','original phase')
%         xlabel('Frequency (Hz)','Fontsize',sn);
%         ylabel('Phase angles','Fontsize',sn);
%         title('compare original phase with phase with noise')
%     subplot(212)
%         stem(f_output, normError,'*')
%         legend('phase error')
%         ylim([0 pi])
%         grid on
%         xlabel('Frequency (Hz)','Fontsize',sn);
%         ylabel('Abs phase error [0 - \pi]','Fontsize',sn);
    
%% (C) ALPC + forward subset selection (SFS)
% Option A: SFS starting from 1 single f
% f1 = f_output(1) %Hz
f1 = 38

TC = -1*ERPt; % see my paper for TC
tmax = 0.15;
[output2] =  Utils.getALPC_SFS1_TC(f1, f_angle, TC, tmax,1);
f_sel = output2.fsel; 
            
        
% use the selected subset
n0 = 8;% the step where MSE starts to ascend
f_subset = f_sel(1: n0+1);
[min_MSE, tau_est, pltOut1]  = Utils.apprent_latency_cycles(f_subset, yy, [], fs, 1, -1*ERPt);
fprintf('sub-system 1...')
min_MSE
tau_est

% (1) use the remaining as subset2
f_subset = f_sel(n0+2: end);
[min_MSE, tau_est, pltOut2]  = Utils.apprent_latency_cycles(f_subset, yy, [], fs, 1, -1*ERPt);
fprintf('sub-system 2...')
min_MSE
tau_est

% (2) select from one f from subset2 
% f2 = 78
% [output2] =  getALPC_SFS1_TC(f2, f_angle, TC, tmax,1);
% f_sel = output2.fsel;                 
% n0 = 3;% the step where MSE starts to ascend
% f_subset = f_sel(1: n0+1);
% [min_MSE, tau_est, pltOut2]  = apprent_latency_cycles(f_subset, yy, [], fs, 1, -1*ERPt);
% fprintf('sub-system 2 ...')
% min_MSE
% tau_est


%% plot two fitting lines (tau) in one figure
x1=pltOut1.x1;
y1=pltOut1.y1;
y1b=pltOut1.y2;
tau_slope1 = pltOut1.tau_slope; 
tt1= pltOut1.tt;
mse1 = pltOut1.mse;
e_sq1 = pltOut1.e_sq;

x2=pltOut2.x1;
y2=pltOut2.y1;
y2b=pltOut2.y2;
tau_slope2 = pltOut2.tau_slope; 
tt2= pltOut2.tt;
mse2 = pltOut2.mse;  
e_sq2 = pltOut2.e_sq;
    

sn = 12;    
figure % (1)
    ax(1)=subplot(211);
        plot(tt1, e_sq1); hold on;
        plot(tt1, mean(e_sq1, 1), '--k','LineWidth',2)
%         legend (num2str(x1))
        legend ([cellstr(num2str(x1));'MPE'])
        ylim([0 2])
        grid on
        xlabel('t (ms)','Fontsize',sn)
        ylabel('PE','Fontsize',sn)
        title (['Estimated subsystem y_1'])
    ax(2)=subplot(212);   
        plot(tt2, e_sq2); hold on;
        plot(tt2, mean(e_sq2, 1), '--k','LineWidth',2)
        legend ([cellstr(num2str(x2));'MPE'])
        ylim([0 2])
        grid on
        xlabel('t (ms)','Fontsize',sn)
        ylabel('PE','Fontsize',sn)
        title (['Estimated subsystem y_2'])
    linkaxes(ax,'x')  

    
figure % (specify the color of each line)
    temp = jet(9+4);% 9+4 Different Colors  
    temp1 = temp([1:9],:);
    temp2 = temp([9+1:end],:);

    ax(1)=subplot(211);
        h1=plot(tt1, e_sq1'); hold on;
        set(h1, {'color'}, num2cell(temp1,2));% 9 Different Colors    
        plot(tt1, mean(e_sq1, 1), '--k','LineWidth',2)
%         legend (num2str(x1))
        legend ([cellstr(num2str(x1));'MPE'])
        ylim([0 2])
        grid on
        xlabel('t (ms)','Fontsize',sn)
        ylabel('PE','Fontsize',sn)
        title (['Estimated subsystem y_1'])
    ax(2)=subplot(212);   
        h2=plot(tt2, e_sq2'); hold on;
        set(h2, {'color'}, num2cell(temp2,2));% 4 Different Colors
        plot(tt2, mean(e_sq2, 1), '--k','LineWidth',2)
        legend ([cellstr(num2str(x2));'MPE'])
        ylim([0 2])
        grid on
        xlabel('t (ms)','Fontsize',sn)
        ylabel('PE','Fontsize',sn)
        title (['Estimated subsystem y_2'])
    linkaxes(ax,'x') 
 
    
 figure
    plot(x1, y1b)
    hold on
    plot(x2, y2b)
    plot(x1, y1,'o')
    plot(x2, y2,'o')
%     legend('1','2','3','4')
%     xlim([0 80]) % Hz
    grid on
    xlabel('Frequency (Hz)','Fontsize',sn);
    ylabel('Phase lag (rad)','Fontsize',sn);
    tao1 = sprintf('%.0f', 1000*tau_slope1);
    tao2 = sprintf('%.0f', 1000*tau_slope2);
    title (['Estimated \tau_1 = ',  tao1, ' ms,', '  \tau_2 = ',tao2, ' ms'])
    
   
 figure
    hold on
    for i=1:length(x1)
        h=plot(x1(i), y1(i),'o','LineWidth',2);
        set(h, {'color'}, num2cell(temp1(i,:),2)); 
    end
    plot(x1, y1b,'k')  
    for i=1:length(x2)
        h=plot(x2(i), y2(i),'o','LineWidth',2);
        set(h, {'color'}, num2cell(temp2(i,:),2)); 
    end
    plot(x2, y2b,'k')
    grid on
    xlabel('Frequency (Hz)','Fontsize',sn);
    ylabel('Phase lag (rad)','Fontsize',sn);
    title (['Estimated \tau_1 = ',  tao1, ' ms,', '  \tau_2 = ',tao2, ' ms'])
%     title (['Estimated \tau_p_1 = ',  tao1, ' ms,', '  \tau_p_2 = ',tao2, ' ms'])   
    
 
    

    