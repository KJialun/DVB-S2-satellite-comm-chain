
addpath(genpath('function'));
%% Parameters
clear all
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]

% phase shift
phase=0;
M = 10; % oversampling factor
fsymb = 2*f_cut; % the symbol frequency 
Tsymb = 1/fsymb; % time between two symbols.
fsampling = M*fsymb; % sampling frequency
Npbsa=[1 2 4 6];

%% CFO
fc = 2e+9;
ppm = fc*1e-6;
df=[0 2 4 6 8 10].*ppm;

for i=1:length(df)
Nbps=2;
dw=2*pi*df(i);
beta = 0.3; % roll-off factor

Nbits = 1000; % bit stream length
bits_tx = randi(2,Nbits,1)-1;


%% Mapping of Binary signal 
if(Nbps > 1)
   modulation = 'qam';
else
   modulation = 'pam';
end
symb_tx = mapping(bits_tx,Nbps,modulation);
%% Upsampling
symb_tx_upsample = upsample(symb_tx,M);
%% Implementation of HRC
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
signal_hrrc_tx_uncoded = conv(symb_tx_upsample,  HRC_time);


%% Noise through the channel uncoded
SNR = 1:0.2:20;
signal_power_uncoded = (trapz(abs(signal_hrrc_tx_uncoded).^2))*(1/fsampling); % total power
Eb = signal_power_uncoded/(Nbits); % bit energy
Eb= Eb/2; % Power of signal is half of its complex envelope's power


for j = 1:length(SNR)
    N0 = Eb/(10.^(SNR(j)/10));
    NoisePower = 2*N0*fsampling;
    noise = sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx_uncoded),1)+1i*randn(length(signal_hrrc_tx_uncoded),1));
    % with noise
    signal_rxt = signal_hrrc_tx_uncoded + noise;

    %without noise
%     signal_rxt = signal_hrrc_tx_uncoded;
%% add CFO
    signal_rx = CFO(signal_rxt,fsampling,dw,phase);
    
    %% matched filter
    signal_hhrc_rx = conv(signal_rx, HRC_time);
    signal_hhrc_rx_trunc = signal_hhrc_rx(taps:end-taps+1);
    %% CFO SHIFT Compensate    
    signal_hhrc_rx_trunc_cfo_shift=Compensate_CFO_shift(signal_hhrc_rx_trunc,fsampling,dw,phase);  

    %% Downsampling
    signal_rx_down = downsample(signal_hhrc_rx_trunc, M);
    signal_rx_down_cfo_compensate = downsample(signal_hhrc_rx_trunc_cfo_shift, M);
    if i==1
        draw_noCFO=signal_rx_down;
    end
    %% Demapping
    if(Nbps~=1)
    uncoded_bits_rx = (demapping(signal_rx_down,Nbps,modulation))';
    uncoded_bits_rx_cfo_compensate = (demapping(signal_rx_down_cfo_compensate,Nbps,modulation))';
    else
    uncoded_bits_rx = (demapping(real(signal_rx_down),Nbps,modulation))'; 
    uncoded_bits_rx_cfo_compensate = (demapping(real(signal_rx_down_cfo_compensate),Nbps,modulation))';   
    end
    
    BER1(j,i) = sum(bits_tx ~= uncoded_bits_rx')/length(uncoded_bits_rx');
    BER2(j,i) = sum(bits_tx ~= uncoded_bits_rx_cfo_compensate')/length(uncoded_bits_rx_cfo_compensate');
end


end

BER1=db(BER1);
BER2=db(BER2);

%% plot CFO shift
figure
subplot(1,3,1)
plot(symb_tx,'o')
title('no cfo constellation (no noise)')
subplot(1,3,2)
plot(draw_noCFO,'o')
title('no cfo constellation (noise)')
subplot(1,3,3)
plot(signal_rx_down,'o')
title('with cfo constellation (noise)')
% subplot(2,2,4)
% plot(signal_rx_down_cfo_compensate,'*')
% title('with cfo shift compensate constellation (noise)')
     

%% plot ber cfo shift compensate
figure
plot(SNR,(BER2(:,1))','-',SNR,(BER2(:,2))','-',SNR,BER2(:,3),'-',SNR,BER2(:,4),'-',SNR,BER2(:,5),'-',SNR,BER2(:,6),'-');
legend('0 ppm','2 ppm','4 ppm','6 ppm', '8 ppm', '10 ppm');
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('BER curves cfo shift compensate')
grid on;
hold on