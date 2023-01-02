%% Projet modulation & coding
% Impact of CFO and CPE on the BER performance
addpath(genpath('function'));

clear; close all;

K = 0.01;
fc = 2e+9;
ppm = fc*1e-6;
CFO_values = [0,50,100]*ppm;
phi0 = 0;

show_ber=0; %% if show_ber =1  , we need change CFO_values to [0]*ppm

if show_ber==1
    iter1=100;
    diff_SNR1=1;
    Nbits = 50000; % bit stream length
else
    iter1=1;
    diff_SNR1=50;
    Nbits = 1000;
end

for iter = iter1:100


f_cut = 1e6/2; % cut off frequency of the nyquist filter [Mhz]
M = 100; % oversampling factor (mettre à 100?)
fsymb = 2*f_cut; % symbol frequency
fsampling = M*fsymb; % sampling frequency
ts = 1/fsampling;
Tsymb = 1/fsymb; % time between two symbols
beta = 0.3; % roll-off factor
Nbps = 1; % number of bits per symbol
if(Nbps > 1)
   modulation = 'qam';
else
   modulation = 'pam';
end
bits_tx = randi(2,Nbits,1)-1;

tshift = 20;

%% Mapping
symbol_tx = mapping(bits_tx,Nbps,modulation);

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
        SNR =1:50;
       for diff_SNR=diff_SNR1:length(SNR)
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
        
        %% Time shift
        symbol_rx_upsampled = signal_rx_CFO(1+tshift:end);
        
        %% Gardner algrithm 
         [corr, error]= GardenerK(symbol_rx_upsampled,K,M);        
         time_error(iter,:,diff_CFO,diff_SNR) = (tshift/M+error(1:end)).'*Tsymb;
         
         k=[0.01,0.02,0.1];
         for i=1:3
         [~, diff_K_error]= GardenerK(symbol_rx_upsampled,k(i),M);   
         diff_K_time_error(iter,:,i,diff_SNR) = (tshift/M+diff_K_error(1:end)).'*Tsymb;
         end
           %% Downsampling
           if show_ber==1
                signal_rx_down = downsample(symbol_rx_upsampled, M);

                %% Demapping
                if(Nbps~=1)
                bits_rx = (demapping(signal_rx_down,Nbps,modulation))';
                bits_rx_COR = (demapping(corr,Nbps,modulation))';
                else
                bits_rx = (demapping(real(signal_rx_down),Nbps,modulation))';
                bits_rx_COR = (demapping(real(corr),Nbps,modulation))';
                end
                BER1(diff_SNR,1) = sum(bits_tx ~= bits_rx')/length(bits_rx');
                BER2(diff_SNR,1) = sum(bits_tx ~= bits_rx_COR')/length(bits_rx_COR');             
           end        
       end
   end
end


%% Plot TIME ERROR results CFO AFFECT
if show_ber==0
SNR=50; % desired SNR=50 db
% load gardnerCFO.mat
time_error_mean = mean(time_error);
time_error_stdv = std(time_error);
mean1 = time_error_mean(1,:,1,SNR);
mean2 = time_error_mean(1,:,2,SNR);
mean3 = time_error_mean(1,:,3,SNR);
stdv1 = time_error_stdv(1,:,1,SNR);
stdv2 = time_error_stdv(1,:,2,SNR);
stdv3 = time_error_stdv(1,:,3,SNR);

figure
plot(mean1,'r');hold on;
plot(mean2,'g')
plot(mean3,'b')
plot(mean1+stdv1,'--r')
plot(mean1-stdv1,'--r')
plot(mean2+stdv2,'--g')
plot(mean2-stdv2,'--g')
plot(mean3+stdv3,'--b')
plot(mean3-stdv3,'--b')

xlabel('Symbols');
ylabel('Time error (mean \pm stdv)');
legend('CFO = 0 ppm','CFO = 50 ppm','CFO = 100 ppm');
title('Convergence of the Gardner algorithm k=0.01')
grid on;

%% Plot Differen K
time_error_mean = mean(diff_K_time_error);
time_error_stdv = std(diff_K_time_error);
mean1 = time_error_mean(1,:,1,SNR);
mean2 = time_error_mean(1,:,2,SNR);
mean3 = time_error_mean(1,:,3,SNR);
stdv1 = time_error_stdv(1,:,1,SNR);
stdv2 = time_error_stdv(1,:,2,SNR);
stdv3 = time_error_stdv(1,:,3,SNR);

figure
plot(mean1,'r');hold on;
plot(mean2,'g')
plot(mean3,'b')
plot(mean1+stdv1,'--r')
plot(mean1-stdv1,'--r')
plot(mean2+stdv2,'--g')
plot(mean2-stdv2,'--g')
plot(mean3+stdv3,'--b')
plot(mean3-stdv3,'--b')

xlabel('Symbols');
ylabel('Time error (mean \pm stdv)');
legend('K = 0.01 ','K= 0.02 ','K = 0.1 ');
title('Convergence of the Gardner algorithm')
grid on;
end

if show_ber==1
    %% after Gardner coorect ber 
figure
plot(SNR,db(BER1(:,1)));
hold on
plot(SNR,db(BER2(:,1)));
legend('without Gardner ','apply Gardner ');
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('BER curves')
end

