%% Projet modulation & coding
% Impact of CFO and CPE on the BER performance
addpath(genpath('function'));
clear; close all;
N_values = [20,30,40];
K_values = [1,8,16];
%% Parameters

for n = 1:length(N_values) 

N = N_values(n);

for iter = 1:500
    
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 100; % oversampling factor (mettre à 100?)
tshift=10;
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
ts = 1/fsampling;
Tsymb = 1/fsymb; % time between two symbols
beta = 0.3; % roll-off factor
Nbps = 2; % number of bits per symbol
modulation = 'qam'; % type of modulation 

Nbits = 100; % data bit stream length
Npilot = N; % Number of pilot bits
bits_pilot = randi(2,1,Npilot)-1;
bits_data = randi(2,1,Nbits)-1;

fc = 2e+9;
ppm = fc*1e-6;
CFO_values = [30]*ppm;
phi0 = 0;

pilot_pos = 30;

%% Mapping of encoded signal
symbol_pilot = mapping(bits_pilot',Nbps,modulation);
symbol_data = mapping(bits_data',Nbps,modulation);
symbol_tx = [symbol_data(1:pilot_pos-1);symbol_pilot;symbol_data(pilot_pos:end)];

%% Upsampling
symbol_tx_upsampled = upsample(symbol_tx,M);

%% Implementation of HHRC
        taps = 365;
        FreqResolution = (1/taps)*fsampling;
        highestfreq = (taps-1)*FreqResolution/2;
        f = linspace(-highestfreq,highestfreq,taps);
        Rc_freq = RC(f,Tsymb,beta); % Raised cosine
        Rc_time = fftshift(ifft(ifftshift(Rc_freq)));
        normal_Rc_time=Rc_time/max(Rc_time);
        normal_Rc_freq=fftshift(abs(fft(normal_Rc_time)));
        HRC_freq = sqrt(normal_Rc_freq);
        HRC_time = fftshift(ifft(ifftshift(HRC_freq)));

%% Convolution
signal_tx = conv(symbol_tx_upsampled, HRC_time);

%% CFO LOOP
   for diff_CFO = 1:length(CFO_values)
%% Noise through the channel
      SNR =1:20;
       for diff_SNR=1:length(SNR)
        signal_power = (trapz(abs(signal_tx).^2))*(1/fsampling); % total power
        Eb = signal_power/(Nbits); % bit energy
        Eb= Eb/2; % Power of signal is half of its complex envelope's power
        N0 = Eb/(10.^(SNR(diff_SNR)/10));
        NoisePower = 2*N0*fsampling;
        noise = sqrt(NoisePower/2)*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));
        signal_rx = signal_tx + noise;
        %% add CFO
        cfo = CFO_values(diff_CFO);
        dw=2*pi*cfo;
        signal_rx_CFO = CFO(signal_rx,fsampling,dw,phi0);
        %% matched filter
        signal_rx_CFO = conv(signal_rx_CFO, HRC_time);
        signal_rx_CFO = signal_rx_CFO(taps:end-taps+1);
        
        signal_rx = conv(signal_rx, HRC_time);
        signal_rx = signal_rx(taps:end-taps+1);
         
        %% Downsampling
        symbol_rx = downsample(signal_rx_CFO(1+tshift:end), M);
%         symbol_rx = [symbol_rx(1:est_n/M);symbol_data(est_n/M:end)];
        NO_cfo_symbol_rx=downsample(signal_rx, M);
        
            for k = 1:length(K_values)
                K = K_values(k);
                %% Frame and frequency acquisition
                [est_n,est_cfo] = cfoEstimate(symbol_rx,symbol_pilot,Tsymb,K);
                [no_cfo_est_n,no_cfo_est_cfo] = cfoEstimate(NO_cfo_symbol_rx,symbol_pilot,Tsymb,K);

        %         exp_cfo2 = exp(-1j*(2*pi*cfo*(0:length(symbol_rx)-1)*ts*M))';
        %         symbol_rx = symbol_rx.*exp_cfo2;

                %% Demapping
                bits_rx = (demapping(symbol_rx,Nbps,modulation))';
        %         BER(j,m) = length(find(bits_data ~= bits_rx))/length(bits_rx');
                time_error(iter,diff_SNR,k,n) = est_n - pilot_pos;
                freq_error(iter,diff_SNR,k,n) = cfo + est_cfo;
                No_cfo_time_error(iter,diff_SNR,k,n) = no_cfo_est_n- pilot_pos;
                No_cfo_freq_error(iter,diff_SNR,k,n) = cfo +no_cfo_est_cfo;
            end
       end
   end
end
end

%% Plot results
% load CFO_K_N.mat
time_error_mean = std(time_error);  
freq_error_mean = std(freq_error);
No_cfo_time_error_mean = std(No_cfo_time_error);  
No_cfo_freq_error_mean = std(No_cfo_freq_error);

N_values = [20, 30, 40];
K_values = [1, 8, 16];
%% no cfo , no time shift  TIME ERROR
figure
plot(SNR,No_cfo_time_error_mean(1,:,1,3),'-r');hold on;
plot(SNR,time_error_mean(1,:,1,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Time error stdv [samples]');
legend('NO_CFO','CFO');
title(' no CFO Time error variances')
grid on;
%% no cfo , no time shift  freq ERROR
figure
plot(SNR,No_cfo_freq_error_mean(1,:,2,1),'-r');hold on;
plot(SNR,freq_error_mean(1,:,2,1),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('FREQ error stdv [samples]');
legend('NO_CFO','CFO');
title(' no CFO FREQ error variances')
grid on;
%% with cfo 
figure
plot(SNR,time_error_mean(1,:,2,1),'-r');hold on;
plot(SNR,time_error_mean(1,:,2,2),'-g');
plot(SNR,time_error_mean(1,:,2,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Time error stdv [samples]');
legend('N = 20, K = 8','N = 30, K = 8','N = 40, K = 8');
title('Time error variances')
grid on;

figure
plot(SNR,time_error_mean(1,:,1,3),'-r');hold on;
plot(SNR,time_error_mean(1,:,2,3),'-g');
plot(SNR,time_error_mean(1,:,3,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Time error stdv [samples]');
legend('N = 30, K = 1','N = 30, K = 8','N = 30, K = 16');
title('Time error variances')
grid on;

figure
plot(SNR,freq_error_mean(1,:,2,1),'-r');hold on;
plot(SNR,freq_error_mean(1,:,2,2),'-g');
plot(SNR,freq_error_mean(1,:,2,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Frequency error stdv [ppm]');
legend('N = 20, K = 8','N = 30, K = 8','N = 40, K = 8');
title('Frequency error variances')
grid on;

figure
plot(SNR,freq_error_mean(1,:,1,3),'-r');hold on;
plot(SNR,freq_error_mean(1,:,2,3),'-g');
plot(SNR,freq_error_mean(1,:,3,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Frequency error stdv [ppm]');
legend('N = 30, K = 1','N = 30, K = 8','N = 30, K = 16');
title('Frequency error variances')
grid on;
