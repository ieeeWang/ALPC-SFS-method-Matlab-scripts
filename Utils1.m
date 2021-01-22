% helper functions for implementing ALPC-SFS
% Lei Wang, ieeewangl@gmail.com
% Sept. 28, 2020 @Radboud Uni.

function Utils1 = Utils1

    Utils1.apprent_latency_cycles = @apprent_latency_cycles;
    Utils1.getfft_compx = @getfft_compx;
    Utils1.getfft_c = @getfft_c;
    
    Utils1.MSPC_2b = @MSPC_2b;
    Utils1.MSPC_latency = @MSPC_latency;
        
    Utils1.getALPC_SFS1_TC = @getALPC_SFS1_TC;
    Utils1.apprent_latency_cycles_angle = @apprent_latency_cycles_angle;
    
    % perform a batch process (ALPC-SFS) for signals with all SNRs
    Utils1.autoSelect_fsubset_f1 = @autoSelect_fsubset_f1;
    
    Utils1.get_angle_outputF = @get_angle_outputF;
end


function [min_mse, tao_est, pltOut] = apprent_latency_cycles(f_output, yy,...
    alfa_in, fs, plt, t_stimuli)
%  estimate time delay - apprent latency from phase coherence.
% OUTPUT
%     min_mse - the minimum MSE
%     tao_est - estimated time delay
% INPUT
%     yy - it MUST be a 12-s time series, the function 'getfft_compx.m' is
%     based on 1/12 f resolustion. Otherwise change this function.
%     tmax - search minimum MES within tmax (s)
%     t_stimuli [s] - at when stimuli starting (e.g., 0.3s)
%     alfa_in - initial phase of output f, it can be [] or a vector (for 
%     phase compensation)

pltOut ={}; % for plot outside this function
tmax = 0.15; % sec
% tmax = 2; % sec 

f_output = f_output(:)';
alfa_0=[];
if isempty(alfa_in)
    alfa_0 = zeros(size(f_output));
    alfa_0(alfa_0<0) = alfa_0(alfa_0<0)+2*pi;
else
    alfa_0 = alfa_in(:)'; % a row
end


% get FFT complex value at 'f_input'
P_compx = [];
for i=1:length(f_output)
    P_compx(i) = getfft_compx(yy, fs, f_output(i)); %N = 33 %Hz
end
% angles lie between [-pi, +pi]
angle_f_out = angle(P_compx);

% map [-pi, +pi] to [0, 2pi]
alfa_f = angle_f_out;
alfa_f(alfa_f<0) = alfa_f(alfa_f<0)+2*pi;

tao = [0: 0.0001 : tmax]; % s <<< search minimum MES within 150 ms
tao = tao + t_stimuli;% s. start from 0.3 s <<<--- for stimuli starting at 0.3s

tt = 1000*tao;% ms
% MSE (rad) for delay y (new!) - exp(1i*rag) (i.e., a vector)
e_sq = [];
for  i=1:length(f_output)    
    tmpL = exp(1i*alfa_f(i));
    tmpR = exp(1i*mod(alfa_0(i)-2*pi*tao*f_output(i), 2*pi));
    e_sq(i,:) = abs(tmpL-tmpR); % within [0 - 2]
end
mse = mean(e_sq, 1);
min_mse = min(mse);
tao_est = tt(mse == min(mse)); % ms
% tao_est = tao_est - 1000*t_stimuli; % <<<<<<<<<<<<<<<<<<<<<<---------
% disp(tao_est)

if plt == 1
    % plot the MSE
    figure
        ax(1)=subplot(211);
            plot(tt, e_sq);
            legend (num2str(f_output'))
            grid on
            ylim([0 2])
            xlabel('t (ms)')
            ylabel('Phase error')
        ax(2)=subplot(212);
            plot(tt, mse);
            grid on
            ylim([0 2])
            xlabel('t (ms)')
            ylabel('MPE')
            linkaxes(ax,'x')
end

    %  plot the slope (f vs. phase lag)
    % (1) get phase lags based on estimate of how many cycles passed
    if isempty(alfa_in) 
        p2 = floor(f_output*tao_est/1000);
        delta_lag = 2*pi*(p2) - alfa_f; % for delay y   
    else
        p2 = floor(f_output*tao_est/1000 + ((alfa_f - alfa_0)/(2*pi)));
        delta_lag = 2*pi*(p2) + alfa_0 - alfa_f; % for delay y (new!)
    end
%     (2) get the 'psudo-accurate' phase lag with 0 errors
%     delta_lag = 2*pi*f_output*tao_est/1000; % 2pi*tao  

    x1 = f_output;
    y1 = delta_lag;
    x1 = x1(:); y1=y1(:);
    order = 1;
    pline = polyfit(x1,y1,order); 
    slope = pline(1);
    tau_slope = slope/(2*pi);
    y2 = polyval(pline, x1);
    
if plt == 1
    % (2) calcute the slope (tau) based on 'unwrapping' points
    figure
        plot(x1, y1,'o')
        hold on 
        plot(x1, y2)
%         xlim([0 100]) % Hz
        grid on
        xlabel('Frequency (Hz)','Fontsize',14);
        ylabel('Phase lag (rad)','Fontsize',14);
        tau_print = sprintf('%.4f', tau_slope);
        title (['\tau (slope/2\pi) = ', tau_print, ' s'])
end       
        
% for plot outside this function    
pltOut.x1=x1;
pltOut.y1=y1;
pltOut.y2=y2;
pltOut.tau_slope = tau_slope; 
pltOut.tt=tt;
pltOut.mse=mse;  
pltOut.e_sq=e_sq;
end 


function [P_compx] = getfft_compx(X, fs, N) 
% get the complex DFT at an integer frequency N, (0<N<fs/2)

% f resolusion =  1/t, t is duration (sec) of X
% DFT is computed without setting a NFFT, see explanation in My_project\test\FFT_with_out_NFFT.m
% OUTPUT
% P_compx - complex value at the integer frequency N (Hz)


Y = fft(X); % same with Y = fft(X, lenth(X)); to ensure the maximum f resolusion
L = length(X); % f resolusion = fs/L = 1/t
Y_comp = Y/L;
P1 = Y_comp(1 : (L/2+1)); % first half including the median point
P1(2:end-1) = 2*P1(2:end-1);
f = fs/2*(0:(2/L):1); % same with f = Fs*(0:(L/2))/L;

% P_compx = P1((f==N)); % it makes error when 99 != 99.0000 
P_compx = P1(12*N+1); % use this when f resolution = 1/12 Hz
end


function [P1, f] = getfft_c(X, fs) 
% get the Amplitude Spectrum of X
% f resolusion =  1/t, t is duration (sec) of X
% DFT is computed without setting a NFFT, see explanation in: 
% My_project\test\FFT_with_out_NFFT.m
% OUTPUT
% P1 - |Y(f)|
% P2 - |Y(f)|^2

Y = fft(X); % same with Y = fft(X, lenth(X)); 
L = length(X);
P_abs = abs(Y/L);
P1 = P_abs(1 : (L/2+1)); % first half including the median point
P1(2:end-1) = 2*P1(2:end-1);
f = fs/2*(0:(2/L):1); % same with f = Fs*(0:(L/2))/L;
% f resolusion = fs/L = 1/t
P2 = P1.^2; 

% figure
%     plot(f, P1) 
%     title('Single-Sided Amplitude Spectrum of X(t)')
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')
end


function [C,C_p,Ang,TimeD,f_sigma]=MSPC_2b(X,Y,fe) 
%X,Y are the fourier transform of the input and the output
%fe is the vector contains fundamental frequency
%f is the whole frequency band
% update by Lei based on 'MSPC_2.m'

k=0;
n = length(fe);
C = zeros(n^2,1);
C_p=zeros(n^2,1);
Ang = zeros(n^2,1);
TimeD = zeros(n^2,1);
f_sigma = zeros(n^2,1);

% fi+fj : 6 elements (righ-up angle)
for i = 1:length(fe)
    for j = i:length(fe)
        k = k+1;
        f_sigma(k)=fe(i)+fe(j)-1;
        CSD_i = X(fe(i),:).*X(fe(j),:).*conj(Y(f_sigma(k),:));
        % normalized in single trials
        CSD_i_normalized = CSD_i./abs(CSD_i);
        C_p(k) = mean (CSD_i_normalized);
        Ang(k) = angle(C_p(k));
        C(k)= abs(C_p(k));
        if Ang(k)<0
            Ang(k)=Ang(k)+2*pi;
        end
        %PhaseD_p_EMGflex_degree_nonlinear(k) = Ang(k)*180/pi;
        %TimeD(k) = (Ang(k)/(2*pi))*1000/(f_sigma(k)-1); % time delay
        %hold on
        %plot(f(f_sigma(k)),C(k),'b*');
    end
end

% fj-fi
for i = 1:length(fe)-1
    for j = (i+1):length(fe)
         k = k+1;
         f_sigma(k)=fe(j)-fe(i)+1;
         CSD_i = X(fe(j),:).*conj(X(fe(i),:)).*conj(Y(f_sigma(k),:));
        
        % normalized in single trials
        CSD_i_normalized = CSD_i./abs(CSD_i);
        C_p(k) = mean (CSD_i_normalized);
        Ang(k)= angle(C_p(k));
        C(k)= abs(C_p(k));
        
        if Ang(k)<0
            Ang(k)=Ang(k)+2*pi;
        end
        %TimeD(k) = (Ang(k)/(2*pi))*1000/(f_sigma(k)-1);
    end
end
% [f_sigma,index]=sort(f_sigma(1:9),'ascend');
[f_sigma,index]=sort(f_sigma,'ascend');
f_sigma = f_sigma(:);

C=C(index);
C_p=C_p(index)./C;
Ang = Ang(index);
TimeD = (Ang./(2*pi))*1000./(f_sigma-1);
end


function [min_MSE, tao_est, pltOut] = MSPC_latency(f_output, C_p_2, plt, mytitle)
%  estimate time delay based on MSPC output and plot the slope line
% INPUT
%     f_output - frequency vec
%     C_p_2 - complex value of the frequency vec
% OUTPUT
%     min_mse - the minimum MSE
%     tao_est - estimated time delay
% <<< search minimum MES within 150 ms
 
pltOut ={};
f_output = f_output(:);
% step = 0.1 ms
% tao = [0: 0.0001 : 0.06]; % s
tao = [0: 0.0001 : 0.15]; % s <<< search minimum MES within 150 ms
tt = 1000*tao;% ms
% MSE (rad) for delay y (new!) - exp(1i*rag) (i.e., a vector)
e_sq = [];
for  i=1:length(f_output)
    e_sq(i,:) = abs(C_p_2(i) - exp(1j*(f_output(i))*2*pi*tao));% within [0 - 2] 
end
mse = mean(e_sq, 1);
min_MSE = min(mse);
tao_est = tt(mse == min(mse)); % ms

%% fit the slope (f vs. phase lag)
% angle2 = atan(imag(C_p_2)./real(C_p_2));
angle2 = angle(C_p_2); % this one seems better than using 'atan'

% correct 2*pi problem: minus phase + 2*pi
angle2(angle2<0) = angle2(angle2<0)+2*pi;

p2 = floor(f_output*tao_est/1000);
p2 =p2(:); angle2=angle2(:);
delta_lag = 2*pi*(p2)  + angle2; % for delay y

x1 = f_output;
y1 = delta_lag;
x1 = x1(:); y1=y1(:);
order = 1;
pline = polyfit(x1,y1,order); 
slope = pline(1);
tau_slope = slope/(2*pi);
y2 = polyval(pline, x1);

%% plot the PE of each freq and latency
% option A: using defualt color
% if plt == 1
%     figure
%     subplot(211);
%         plot(tt, e_sq); hold on;
%         plot(tt, mse,'-.k','LineWidth',2);
%         %legend (num2str(f_output))
%         legend ([cellstr(num2str(f_output));'MPE'])
%         grid on
%         ylim([0 2])
%         xlabel('t (ms)','Fontsize',12)
%         ylabel('Phase error (PE)','Fontsize',12)
%         title(mytitle)
%     % (2) calcute the slope (tau) based on 'unwrapping' points
%     subplot(212)
%         plot(x1, y1,'o')
%         hold on
%         plot(x1, y2)
% %         xlim([0 100]) % Hz
%         grid on
%         xlabel('Frequency (Hz)','Fontsize',12);
%         ylabel('Phase lag (rad)','Fontsize',12);
%         %title (['MSPC: \tau (slope/2\pi) = ', num2str(tau_slope) , ' s'])
%         title (['\tau (d_\phi/2\pi) = ', sprintf('%.3f',tau_slope), ' s'])
% end        
 
% option B: using jet color
if plt == 1
    n_freq = size(e_sq,1);
    temp_color = jet(n_freq);% 9+4 Different Colors  
    figure
    subplot(211);
        h1= plot(tt, e_sq); hold on;
        set(h1, {'color'}, num2cell(temp_color,2));% Different Colors
        plot(tt, mse,'-.k','LineWidth',2);
        %legend (num2str(f_output))
        legend ([cellstr(num2str(f_output));'MPE'],'NumColumns',2)
        grid on
        ylim([0 2])
        xlabel('t (ms)','Fontsize',12)
        ylabel('Phase error (PE)','Fontsize',12)
        title(mytitle)
    % (2) calcute the slope (tau) based on 'unwrapping' points
    subplot(212)
        plot(x1, y1,'o')
        hold on
        plot(x1, y2)
%         xlim([0 100]) % Hz
        grid on
        xlabel('Frequency (Hz)','Fontsize',12);
        ylabel('Phase lag (rad)','Fontsize',12);
        %title (['MSPC: \tau (slope/2\pi) = ', num2str(tau_slope) , ' s'])
        title (['\tau (d_\phi/2\pi) = ', sprintf('%.1f',1000*tau_slope), ' ms'])
end 

% sprintf('%.3f',tau_slope)
% sprintf('%.2f%%',100*recall2)

% for plot [freq vs latency] outside this function    
pltOut.x1=x1;
pltOut.y1=y1;
pltOut.y2=y2;
pltOut.tau_slope = tau_slope; 
pltOut.tt=tt;
pltOut.mse=mse;   
pltOut.e_sq=e_sq;
end


function [output2] =  getALPC_SFS1_TC(f1, f_angle, TC, tmax, plt)
% ALPC-based SFS starting from 1 single f, ALPC with time composation (TC)
% INPUT
% TC: time compensation, e.g.,
%     TC = -0.3s in this EEG study
%     TC = 0s in the simulation
%     f_angle - [f_whole, phase_whole]
%     tmax [s] - e.g., tmax = 0.15 or 0.05 sec, [0 tmax] is the potential
%     range of |tau|.
% OUTPUT
%     output2 
% Lei@RU Sept19 

% (1) initial values 
f_whole = f_angle(:,1);
phase_whole = f_angle(:,2);
[fset, ~] = getRemainingVec2(f_whole, f1);

f_sel = f1;
f_left = fset;

% (2) forward
stepR = {}; min_mse_step = []; tau_step =[];
for i = 1:length(fset)
    tmp_mse =[];
    N= length(f_left);
    for j = 1:N
        f_test = [f_sel, f_left(j)]; % new candidate set: 'f_sel'
        P_angle = get2ndcolumn_item(f_angle, f_test);
        [tmp_mse(j,1), tmp_mse(j,2), ~] = apprent_latency_cycles_angle(...
            f_test, P_angle, [], 0, TC, tmax);
    end
    stepR{i} = tmp_mse;
    
    [output]= select_nextF2(f_left, tmp_mse, TC);
    f_sel = [f_sel, output.fsel]; % add selected f in selected set
    f_left(f_left==output.fsel)=[];% remove selected f in left set
    min_mse_step(i) = output.min_mse;% minimun MPE on each step
    tau_step(i) = output.tau_sel;% estimated tao on each step
end

output2 ={};
output2.fsel = f_sel;%the order of f baased on SFS
output2.stepR = stepR;
output2.MPE = min_mse_step;
output2.tau_step = tau_step; % corresponding tau of the one with min MPE

% visualize SFS 
if plt
    sn = 12;
    figure
        for i=1:length(stepR)
            temp = stepR{i};
            temp = temp(:,1); 
            subplot(211);
            plot(i*ones(size(temp)), temp, 'ob')
            hold on
        end
        plot(min_mse_step,'--r')
        grid on
            xlabel('Step','Fontsize',sn)
            ylabel('MPE','Fontsize',sn)
            %title(['SFS starting from ',num2str(f1),' Hz'])
            title(['SFS on each step'])
        for i=1:length(stepR)
            temp = stepR{i};
            temp = temp(:,2); % latencies
            subplot(212);
            plot(i*ones(size(temp)), temp-TC*1000, 'ob') %TC: +0.3 s in plot
            hold on
        end
        %plot(tau_step,'--r')
        plot(tau_step-TC*1000,'--r')%TC: +0.3 s in plot
        grid on
    %     ylim([0 100])
            xlabel('Step','Fontsize',sn)
            ylabel('\tau (ms)','Fontsize',sn)
end
end


function [fset, loc] = getRemainingVec2(f_all, f1)
% loc is the location index of f1 in f_all

loc = [];
k = 0;
for i=1:length(f_all)
    if  sum(f_all(i)==f1)>0
        k=k+1;
        loc(k) = i;
    end
end
f_all(loc)=[];
fset = f_all;
end


function P_angle = get2ndcolumn_item(f_angle, f_test)
% f_angle - two-column items
% f_test - items in 1st column
% P_angle - corresponding items in 2nd column 

c1 = f_angle(:,1);
c2 = f_angle(:,2);
loc = [];
for i = 1:length(f_test)
    loc(i) = find(c1 == f_test(i));
end
P_angle = c2(loc);
end


function [output] = select_nextF2(f_left, MPE_tau, TC)
% Here set the rules of SFS to select the next f and update 'fsel'
% without using a termination rule for SFS
% MPE_tau 
%     - [min_MSE, tau_est]
% f_left
%     - A vector of candidate f
% TC: time compensation, e.g.,
%     TC = -0.3s in this EEG study
%     TC = 0s in the simulation
% Lei@RU Sept19  

R = [3, 150]; % ms, the resonable range of 3<tau<150 ms
R = R + TC*1000;% ms, for time compensation

% (1) determine a subset whose tau is within the predefined range: R
mse = MPE_tau(:,1);
tau = MPE_tau(:,2);
loc = (tau>R(1)) & (tau<R(2)); 

% update according the subset
f2 = f_left(loc);
mse2 = mse(loc);
tau2 = tau(loc);

% (2) select the f with the minimum MSE
min_mse = min(mse2);
loc2 = (mse2 == min_mse);
fsel = f2(loc2);
tau_sel = tau2(loc2);
if length(fsel)>1
    fsel =fsel(1); % chose only 1st one
end
if length(tau_sel)>1
    tau_sel =tau_sel(1); % chose only 1st one
end

% (3) output
output = {};
output.fsel = fsel;% the selected f
output.min_mse = min_mse; % mse of the selected f
output.tau_sel = tau_sel; % tau of the selected f

end



function [min_mse, tao_est, pltOut] = apprent_latency_cycles_angle(f_output, ...
    angle_f_out, alfa_in, plt, t_compens, tmax)
%  estimate time delay - apprent latency based on preceding cycles or phase
%  coherence.
% OUTPUT
%     min_mse - the minimum MSE
%     tao_est - estimated time delay
% INPUT
%    P_compx - complex vector based on FFT on f_output.
%     tmax - search minimum MES within tmax (s)
%     t_compens [s] - time conpensation for aligning stimuli and EEG (e.g.,
%     -0.3s), see more details in my paper draft.
%     alfa_in - initial phase of output f, it can be [] or a vector
% to directly use normorlized phase (P_compx), e.g., AVG of 1-s epochs

pltOut ={}; % for plot outside this function
% tmax = 0.15; % sec
% tmax = 2; % sec 

f_output = f_output(:)';
alfa_0=[];
if isempty(alfa_in)
    alfa_0 = zeros(size(f_output));
    alfa_0(alfa_0<0) = alfa_0(alfa_0<0)+2*pi;
else
    alfa_0 = alfa_in(:)'; % a row
end


% (option A.1) get FFT complex value at 'f_input'
% P_compx = [];
% for i=1:length(f_output)
%     P_compx(i) = getfft_compx(yy, fs, f_output(i)); %N = 33 %Hz
% end

% (option A.2) get normorlized phase of 12 epochs (1s)
% output = getPC_12_1s_epoch(yy, fs, f_output);
% P_compx = output.Phase_target;

% angles lie between [-pi, +pi]
% angle_f_out = angle(P_compx);
angle_f_out = angle_f_out(:)';

% map [-pi, +pi] to [0, 2pi]
alfa_f = angle_f_out;
alfa_f(alfa_f<0) = alfa_f(alfa_f<0)+2*pi;

tao = [0: 0.0001 : tmax]; % s <<< search minimum MES within 150 ms
% tao = [0: 0.0001 : tmax(2)]; % tmax in SEC

tao = tao + t_compens;% s

tt = 1000*tao;% ms
% MSE (rad) for delay y (new!) - exp(1i*rag) (i.e., a vector)
e_sq = [];
for  i=1:length(f_output)    
    tmpL = exp(1i*alfa_f(i));
    tmpR = exp(1i*mod(alfa_0(i)-2*pi*tao*f_output(i), 2*pi));
    e_sq(i,:) = abs(tmpL-tmpR); % within [0 - 2]
end
mse = mean(e_sq, 1);
min_mse = min(mse);
tao_est = tt(mse == min_mse); % ms

tao_est = tao_est(1); % if more than one minimum exist, use 1st

if plt == 1
    %% plot the MSE
    figure
    sn=12;
        ax(1)=subplot(211);
            plot(tt, e_sq);
            legend (num2str(f_output'))
            grid on
            ylim([0 2])
            xlabel('t (ms)','Fontsize',sn)
            ylabel('PE','Fontsize',sn)
            title('ALPC with SFS')
        ax(2)=subplot(212);
            plot(tt, mse);
            grid on
            ylim([0 2])
            xlabel('t (ms)','Fontsize',sn)
            ylabel('MPE','Fontsize',sn)
            linkaxes(ax,'x')
          
end

    %%  plot the slope (f vs. phase lag) 
    % (1) get phase lags based on estimate of how many cycles passed
    if isempty(alfa_in) 
        p2 = floor(f_output*tao_est/1000);
        delta_lag = 2*pi*(p2) - alfa_f; % for delay y   
    else
        p2 = floor(f_output*tao_est/1000 + ((alfa_f - alfa_0)/(2*pi)));
        delta_lag = 2*pi*(p2) + alfa_0 - alfa_f; % for delay y (new!)
    end
    % (2) get the 'psudo-accurate' phase lag with 0 errors
%     delta_lag = 2*pi*f_output*tao_est/1000; % 2pi*tao      
    

    x1 = f_output;
    y1 = delta_lag;
    x1 = x1(:); y1=y1(:);
    order = 1;
    pline = polyfit(x1,y1,order); 
    slope = pline(1);
    tau_slope = slope/(2*pi);
    y2 = polyval(pline, x1);
    
if plt == 1
    % (2) calcute the slope (tau) based on 'unwrapping' points
    figure
        plot(x1, y1,'o')
        hold on 
        plot(x1, y2)
%         xlim([0 100]) % Hz
        grid on
        xlabel('Frequency (Hz)');% ,'Fontsize',14
        ylabel('Phase lag (rad)');
        title (['\tau (slope/2\pi) = ', num2str(tau_slope) , ' s'])
end       
        
% for plot outside this function    
pltOut.x1=x1;
pltOut.y1=y1;
pltOut.y2=y2;
pltOut.tau_slope = tau_slope; 
pltOut.tt=tt;
pltOut.mse=mse;  
pltOut.e_sq=e_sq;
    
end


function output = autoSelect_fsubset_f1(f1, f_output, angle_Y0)
% perform a batch process (ALPC-SFS) for signals with all SNRs

output ={};
TC = 0; % see my paper for TC
tmax = 0.060; % sec
% f1 = 10; % Hz % f_output = [4;6;10;34;38;42;44;48;54]
% f1 = 34; % Hz % f_output = [4;6;10;34;38;42;44;48;54]

plt=0;
f_select ={}; MPE_sel = []; tau_sel = [];
for j=1:size(angle_Y0, 2)
  
    angle_temp = angle_Y0(:,j);

    % SFS on all steps
    [min_MSE, tau_est, pltOut1] = apprent_latency_cycles_angle(...
                f_output, angle_temp, [], plt, TC, tmax); 

    f_phase = [f_output, angle_temp(:)];
    [output2] =  getALPC_SFS1_TC(f1, f_phase, TC, tmax,plt);

    f_sel = output2.fsel;
    MPE_step = output2.MPE;
    tau_step = output2.tau_step;

    % apply the termination rule: (1) dMPE <0.1 (2) MPE<0.1    
    n0=1;
    for i=2:length(MPE_step)
        dMPE = MPE_step(i)-MPE_step(i-1);
        MPE = MPE_step(i);
        if MPE <0.2 && dMPE<0.1
            n0=i;
        else
            break
        end
    end
    % use only selected f subset
    f_subset = f_sel(1: n0+1);
    P_angle = get2ndcolumn_item(f_phase, f_subset);
    [min_MSE, tau_est, pltOut1] = apprent_latency_cycles_angle(...
                f_subset, P_angle, [], plt, TC, tmax);
    f_select{j,1} = f_subset;
    MPE_sel(j)= min_MSE;
    tau_sel(j)= tau_est;
end

output.f_select =f_select;
output.MPE_sel =MPE_sel;
output.tau_sel =tau_sel;
end



function angle_Y0 = get_angle_outputF(f_output, YY, fs)
% INPUT
%   f_output - a vec of frequency
% YY - each row is 12-s signal
% OUTPUT
%   angle_Y0 - each column is angle of corresponding frequency

angle_Y0 = []; % angle of all sum-up Y0 with diff 'Pr'
for j=1:size(YY,1)
    yy = YY(j,:); 
    compx_temp=[];
    for i=1:length(f_output)
        compx_temp(i) = getfft_compx(yy, fs, f_output(i)); 
    end
    angle_Y0(:,j) = angle(compx_temp); % angles lie between [-pi, +pi]
end
end

