% %% Projet modulation & coding
addpath(genpath('function'));
addpath(genpath('../'));
%% Parameters
clear 
close
clc
f_cut = 1e6/2; % cut off frequency of the nyquist filter [Mhz]

% phase shift
M = 100; % oversampling factor
fsymb = 2*f_cut; % the symbol frequency 
fsampling = M*fsymb; % sampling frequency
Tsymb = 1/fsymb; % time between two symbols.
Npbsa=[1 2 4 6];
shifts=[0,2,10,20];
for i=1:length(shifts)
        disp(i)
        Nbps=2;
        shift=shifts(i);
        dw=0;
        beta = 0.3; % roll-off factor
        
        Nbits = 1200000; % bit stream length
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
        signal_hrrc_tx_uncoded = conv(symb_tx_upsample, HRC_time);
        
        
        %% Noise through the channel uncoded
        SNR = 1:1:20;
        signal_power_uncoded = (trapz(abs(signal_hrrc_tx_uncoded).^2))*(1/fsampling); % total power
        Eb = signal_power_uncoded/(Nbits); % bit energy
        Eb= Eb/2; % Power of signal is half of its complex envelope's power
        
        
        for j = 1:length(SNR)
            N0 = Eb/(10.^(SNR(j)/10));
            NoisePower = 2*N0*fsampling;
            noise = sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx_uncoded),1)+1i*randn(length(signal_hrrc_tx_uncoded),1));
        
            signal_rxt = signal_hrrc_tx_uncoded + noise;
            %signal_rx = CFO(signal_rxt,fsampling,dw,0);
        
            %% matched filter
            signal_hhrc_rx = conv(signal_rxt, HRC_time);
            signal_hhrc_rx_trunc = signal_hhrc_rx(taps:end-taps+1);
            %% shift added
            signal_hhrc_rx_trunc=signal_hhrc_rx_trunc(1+shift:end,:);

            %% Downsampling

            signal_rx_down = downsample(signal_hhrc_rx_trunc, M);

            %% Demapping
            if(Nbps~=1)
            bits_rx = (demapping(signal_rx_down,Nbps,modulation))';
            else
            bits_rx = (demapping(real(signal_rx_down),Nbps,modulation))';
            end
            BER1(j,i) = sum(bits_tx ~= bits_rx')/length(bits_rx');
        end
        scatterData(:,i) = signal_rx_down;
end
BER1=20*log(BER1);
%% Plot BER results
figure
for i=1:width(BER1)
plot(SNR,(BER1(:,i)));
hold on
end
legend('no time shift','t0=0.02T','t0=0.1T','t0=0.2T', '4', '5','6','7','8','9');
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('BER curves')
grid on;

%% Plot Constellation results for SNR = 20 
scatterplot(scatterData(:,1),1,0,'r.')          
title('Time Shift t_0 = 0T')
grid on
 
scatterplot(scatterData(:,2),1,0,'r.')     
title('Time Shift t_0 = 0.02T')
grid on
    
scatterplot(scatterData(:,3),1,0,'r.')         
title('Time Shift t_0 = 0.1T')
grid on
      
scatterplot(scatterData(:,4),1,0,'r.')         
title('Time Shift t_0 = 0.2T')
grid on
