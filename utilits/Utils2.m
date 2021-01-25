% helper functions for implementing ALPC-SFS
% Lei Wang, ieeewangl@gmail.com
% Sept. 28, 2020 @Radboud Uni.

function Utils2 = Utils2
    Utils2.plot2fitlines_inOne = @plot2fitlines_inOne;
    Utils2.plot_MSPC_ALPC_inOne= @plot_MSPC_ALPC_inOne;
    Utils2.get_2o_3o_output_f = @get_2o_3o_output_f;
    Utils2.get_theory_2nd_order_outputs = @get_theory_2nd_order_outputs;
end


function [output_2nd, output_3rd]= get_2o_3o_output_f(input)
%% by theory
output_2nd = get_theory_2nd_order_outputs(input);
% disp('2nd order outputs:')
% output_2nd'

output_3rd = get_theory_3rd_order_outputs(input);
% disp('3rd order outputs:')
% output_3rd'
end


function f_sigma = get_theory_2nd_order_outputs(input)
% estimating outputs of 2nd order, given input frequenceis

fe = input+1;
k=0;
order=2;
f_sigma = [];

% fi+fj : 6 elements (righ-up angle)
for i = 1:length(fe)
    for j = i:length(fe)
        k = k+1;
        f_sigma(k)=fe(i)+fe(j)-1;
    end
end

% fj-fi
for i = 1:length(fe)-1
    for j = (i+1):length(fe)
         k = k+1;
         f_sigma(k)=fe(j)-fe(i)+1;
    end
end

% ensure positive freq
f_sigma = abs(f_sigma);
[f_sigma, ~]=sort(f_sigma,'ascend');
f_sigma = f_sigma(:)-1;

end


function f_sigma = get_theory_3rd_order_outputs(input)
% estimating outputs of 2nd order, given input frequenceis

fe = input+1; % +1 is to be idexable to its fft x-axis (with a DC)
k=0;
order=3;
f_sigma = [];

%% 2fi+fj: 9=[3*3] elements (including diagonal elements 3*fi)
for i = 1:length(fe)
    for j = 1:length(fe)
        k = k+1;
        f_sigma(k)=2*fe(i)+fe(j)-2;
    end
end
        
%% |2fi-fj|: 9=[3*3] elements (including diagonal elements 2fi-fi=fi)
for i = 1:length(fe)
%     for j = (i+1):length(fe)
    for j = 1:length(fe)
         k = k+1;
         f_sigma(k)= 2*fe(i)-fe(j);
         if f_sigma(k)-1<0
             f_sigma(k)= -(f_sigma(k)-1)+1;
         end
    end
end         
         
%% choose 3 without repeating from input frequencies
Cob3 = nchoosek(fe,3);
for i=1:size(Cob3,1)  
    temp = Cob3(i,:);
    k = k+1; % f1+f2+f3
    f_sigma(k)=temp(1)+temp(2)+temp(3)-2; 
end                                   

Cob3 = nchoosek(fe,3);
for i=1:size(Cob3,1)  
    temp = Cob3(i,:);
    k = k+1; % f1+f2-f3
    f_sigma(k)=temp(1)+temp(2)-temp(3); 
    k = k+1; % f1+f3-f2
    f_sigma(k)=temp(1)+temp(3)-temp(2); 
    k = k+1; % f2+f3-f1
    f_sigma(k)=temp(3)+temp(2)-temp(1); 
end 


% ensure positive freq
f_sigma = abs(f_sigma);
% returns the same data as in A, but with no repetitions. 
f_sigma = unique(f_sigma); %C is in sorted order.
% [f_sigma, ~]=sort(f_sigma,'ascend');

f_sigma = f_sigma(:)-1;
end


function plot2fitlines_inOne(pltOut1, pltOut2)
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

        
figure
    plot(x1, y1b)
    hold on
    plot(x2, y2b)
    plot(x1, y1,'o')
    plot(x2, y2,'o')
    legend('\tau_1','\tau_2','freq in y_1','freq in y_2')
%     xlim([0 80]) % Hz
    grid on
    xlabel('Frequency (Hz)','Fontsize',sn);
    ylabel('Phase lag (rad)','Fontsize',sn);
    tau_print1 = sprintf('%.1f', 1000*tau_slope1);
    tau_print2 = sprintf('%.1f', 1000*tau_slope2);
    title (['Estimated \tau_1 = ', tau_print1, ' ms,', '  \tau_2 = ', tau_print2, ' ms'])
    
end    

function plot_MSPC_ALPC_inOne(pltOut1_alpc, pltOut2_alpc, pltOut1_mspc, pltOut2_mspc)
%% ALPC
x1=pltOut1_alpc.x1;
y1=pltOut1_alpc.y1;
y1b=pltOut1_alpc.y2;
tau_slope1 = pltOut1_alpc.tau_slope; 
tt1= pltOut1_alpc.tt;
mse1 = pltOut1_alpc.mse;
e_sq1 = pltOut1_alpc.e_sq;

x2=pltOut2_alpc.x1;
y2=pltOut2_alpc.y1;
y2b=pltOut2_alpc.y2;
tau_slope2 = pltOut2_alpc.tau_slope; 
tt2= pltOut2_alpc.tt;
mse2 = pltOut2_alpc.mse;  
e_sq2 = pltOut2_alpc.e_sq;

figure
sn=12;
subplot(211)
    plot(x1, y1b)
    hold on
    plot(x2, y2b)
    plot(x1, y1,'o')
    plot(x2, y2,'o')
    grid on
    xlabel('Frequency (Hz)','Fontsize',sn);
    ylabel('Phase lag (rad)','Fontsize',sn);
    tau_print1 = sprintf('%.1f', 1000*tau_slope1);
    tau_print2 = sprintf('%.1f', 1000*tau_slope2);
    textA = ['\tau_1 = ', tau_print1, ' ms'];
    textB = ['\tau_2 = ', tau_print2, ' ms'];
    legend({textA, textB,'freq selected in S_1','freq selected in S_2'},...
        'Location','northwest')
    title (['ALPC-SFS'],'Fontsize',sn)
    

%% MSPC
x1=pltOut1_mspc.x1;
y1=pltOut1_mspc.y1;
y1b=pltOut1_mspc.y2;
tau_slope1 = pltOut1_mspc.tau_slope; 
tt1= pltOut1_mspc.tt;
mse1 = pltOut1_mspc.mse;
e_sq1 = pltOut1_mspc.e_sq;

x2=pltOut2_mspc.x1;
y2=pltOut2_mspc.y1;
y2b=pltOut2_mspc.y2;
tau_slope2 = pltOut2_mspc.tau_slope; 
tt2= pltOut2_mspc.tt;
mse2 = pltOut2_mspc.mse;  
e_sq2 = pltOut2_mspc.e_sq;

subplot(212)
    plot(x1, y1b)
    hold on
    plot(x2, y2b)
    plot(x1, y1,'o')
    plot(x2, y2,'o')
    legend({'\tau_1','\tau_2','freq in y_1','freq in y_2'},'Location','northwest')
    grid on
    xlabel('Frequency (Hz)','Fontsize',sn);
    ylabel('Phase lag (rad)','Fontsize',sn);
    tau_print1 = sprintf('%.1f', 1000*tau_slope1);
    tau_print2 = sprintf('%.1f', 1000*tau_slope2);
    textA = ['\tau_1 = ', tau_print1, ' ms'];
    textB = ['\tau_2 = ', tau_print2, ' ms'];
    legend({textA, textB,'2^{nd}-order freq','3^{rd}-order freq'},...
        'Location','northwest')
    title (['MSPC'],'Fontsize',sn)
    
end