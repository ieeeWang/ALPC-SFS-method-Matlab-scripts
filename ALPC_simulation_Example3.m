clc; clear; close all;
% for paper simulation v3
% generate 2-order responses with different tatencies  
% estimating two tatencies by using mixed phases
% get results of different amplitude ratios
% add noise on Y

% declare the helper function
utils_1 = Utils1; 
utils_2 = Utils2; 
utils_3 = Utils3; 

%% generate signals
fs = 2000;
phase0 = 0.5*pi; %
ERPt = 0; % for simulation

%% Input & output
order = 2;
w1 = [17, 21, 27]; %Hz
tau1 = 0.020; % [sec]
tau2 = 0.015; % [sec]

plt = 0;
% subsys 1: tau1
[x1, y1, t] = generateSignal_12sec(w1, order, tau1, phase0, [], fs, ERPt);
% subsys 2: tau2
[x2, y2, ~] = generateSignal_12sec(w1, order, tau2, phase0, [], fs, ERPt);

% Pr = power ratio of Y1/Y2
Pr = [1/10000, 1/100, 1/36, 1/16, 1/9, 1/4, 1/2,  1,  2, 4, 9, 16, 36, 100, 10000]';
Ar = sqrt(Pr);% amplitude ratio of Y1/Y2
X = x1; % X1 = X2 i.e., two subsys have the same input
Y0 = Ar*y1 + y2; % zero noise

% add noise
SNR_Y = [5, 0, -5, -10, -15, -20]; % dB
% Y with noise
Yn_set={};
for j=1:length(SNR_Y)
    Yn=[]; 
    for i=1:size(Y0,1)
        Yn(i,:) = awgn1(Y0(i,:), SNR_Y(j), 'measured','db');
    end
    Yn_set{j}=Yn;
end

%% phase values on each system output
f_output = utils_2.get_theory_2nd_order_outputs(w1);

% (1) each output with Pr
YY = Y0;
angle_Y0 = utils_1.get_angle_outputF(f_output, YY, fs);
 
% (1) each output (with noise) with Pr
angle_Yn_set={};
for i=1:length(Yn_set)
    YY = Yn_set{i};
    angle_Yn_set{i} = utils_1.get_angle_outputF(f_output, YY, fs);
end
 

% (2) each output in original signals
Y12 = [y1; y2]; % original signals
% angle of all sum-up Y0 with diff 'Pr'
angle_Y12 = utils_1.get_angle_outputF(f_output, Y12, fs);
      
    
%% plot their phases
figure
sn=12;
    stem(f_output, angle_Y0)
    %ylim([-pi pi])
    grid on
    legend(num2str(Ar))
    hold on
    plot(f_output, angle_Y12, '--*')
    hold off
    xlabel('Frequency (Hz)','Fontsize',sn);
    ylabel('Phase angle (rad)','Fontsize',sn);
    title('Output frequency phases on each amplitude ratio')

%% ALPC-SFS
f1 = 10; % Hz % f_output = [4;6;10;34;38;42;44;48;54]
% perform a batch process for signals with all SNRs
output1 = utils_1.autoSelect_fsubset_f1(f1, f_output, angle_Y0);


fprintf('the selected sub-system ...\n') 
f_1 = output1.f_select
MPE_1 = output1.MPE_sel
tau_1 = output1.tau_sel

MPE_noise=[];tau_noise=[];
for i=1:length(angle_Yn_set)
    output_temp = utils_1.autoSelect_fsubset_f1(f1, f_output, angle_Yn_set{i});
    MPE_noise(i,:) = output_temp.MPE_sel; 
    tau_noise(i,:) = output_temp.tau_sel;
end


%% estimated tau
xx = log(Ar);
[fitresult1, gof] = utils_3.createFit_sigmoid_tau(xx(:), tau_1(:));

figure
subplot(211) 
    plot(xx, tau_1, 'o')
    hold on
    plot(xx, tau_noise', '*')
%     plot(fitresult1, xx, tau_1);
    plot(fitresult1);
    hold off
    grid on
%     ylim([20, 52])
%     ylim([15, 20])
    legend('pure signal','5 dB','0 dB','-5 dB','-10 dB','-15 dB','-20 dB','sigmoid fit line')
    xlabel('log(\xi)','Fontsize',sn);
    ylabel('\tau (ms)','Fontsize',sn);
    title('Estimated latencies on mixed signals')
subplot(212) 
    stem(xx, MPE_1, 'o')
    hold on
    stem(xx, MPE_noise', '*')
    legend('pure signal','5 dB','0 dB','-5 dB','-10 dB','-15 dB','-20 dB')
    hold off
    grid on
    xlabel('log(\xi)','Fontsize',sn);
    ylabel('MPE','Fontsize',sn);

%% re-plot above figure   
figure
% pointype = {'.','+','*','s','^','d'}; % 6 SNRs
pointype = {'+','s','.','^','*','d'}; % 6 SNRs
subplot(211) 
    %plot(xx, tau_1, 'o', 'LineWidth',2)
    plot(xx, tau_1, 'o')
    hold on
    for i=1:length(pointype)
        plot(xx, tau_noise(i,:), pointype{i});
    end
    plot(fitresult1);
    hold off
    grid on
%     ylim([15, 20])
    legend('no noise','5 dB','0 dB','-5 dB','-10 dB','-15 dB','-20 dB','sigmoid fitting')
    xlabel('log(\xi)','Fontsize',sn);
    ylabel('\tau (ms)','Fontsize',sn);
    title('Estimated latencies on mixed signals')
subplot(212) 
    plot(xx, MPE_1, '--o')
    hold on
    %plot(xx, MPE_noise', '--*')
    for i=1:length(pointype)
        plot(xx, MPE_noise(i,:), ['--',pointype{i}]);
    end
    legend('no noise','5 dB','0 dB','-5 dB','-10 dB','-15 dB','-20 dB')
    hold off
    grid on
    xlabel('log(\xi)','Fontsize',sn);
    ylabel('MPE','Fontsize',sn);

    
%% sigmoid curve fitting on the normalized tau
% (1) get norm (X, Y) 
X = log(Ar);
yy_tau = tau_1; % or tau_2
L= 1000*tau2; H= 1000*tau1; 
Y = (yy_tau-L)/(H-L);
Y = Y(:);

% (2) get the curve fitting
[fitresult, gof] = utils_3.createFit_sigmoid_norm(X, Y)


return
% function of the nonlinear coeef: y=1./(1+x.^b)
b=-1.17;
x= 0.01*([1:1:1000]);
y=1./(1+x.^b);
figure
    plot(x,y)
    xlabel('r','Fontsize',sn);
    ylabel('coef','Fontsize',sn);    








