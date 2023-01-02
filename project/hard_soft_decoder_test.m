% %% Projet modulation & coding
addpath(genpath('function'));

%% Parameters
clear all

decoding=1; %1=hard 0=soft  the switch of hard decoding and soft decoding 

f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 10; % oversampling factor
fsymb = 2*f_cut; % the symbol frequency 
fsampling = M*fsymb; % sampling frequency
Tsymb = 1/fsymb; % time between two symbols.
Npbsa=[1 2 4 6];
beta = 0.3; % roll-off factor
Npackets = 100; % choose suitable value make sure that Nbits/Nbps is an integer
packetLength = 128;
codedWordLength = 2*packetLength; % make sure rate=1/2  bits_tex/codeword = 1/2
Nbits = Npackets*packetLength; % bit stream length
NcodedBits = Npackets*codedWordLength; % full coded word length
rate=1/2;
bits_tx = randi(2,Nbits,1)-1; % generate random bit stream 0101110011,,,,,, 
bits_tx_coded = zeros(NcodedBits,1);
%% main loop
    Nbps=1;  %%  the simplified algorithm of soft decoding ,  only work for BPSK , we should choose NBPS=1
    %%          for hard decoding ,we could test different Nbps
    %% chanl coding      
    H0 = makeLdpc(packetLength, codedWordLength, 0, 1, 3);
for k=1:Npackets
    packet_tx = bits_tx(1+(k-1)*packetLength : k*packetLength);% choose each packet bits_tex value , each packet length=128
    [codedbits, H] = makeParityChk(packet_tx , H0, 0);  %% encoding 
    bits_tx_coded(1+(k-1)*codedWordLength : k*codedWordLength) = [codedbits packet_tx];
end  
    %% Mapping of Binary signal 
    if(Nbps > 1)
       modulation = 'qam';
    else
       modulation = 'pam';
    end
    symb_tx = mapping(bits_tx_coded,Nbps,modulation);
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
    signal_hrrc_tx_coded = conv(symb_tx_upsample,  HRC_time);
    %% Noise through the channel uncoded
    SNR = 1:1:20;
    BER1 = zeros(length(SNR),5);
    signal_power_coded = (trapz(abs(signal_hrrc_tx_coded).^2))*(1/fsampling); % total power
    Eb = signal_power_coded/(Nbits); % bit energy
    Eb= Eb/2; % Power of signal is half of its complex envelope's power


        for j = 1:length(SNR)
            N0 = Eb/(10.^(SNR(j)/10));
            NoisePower = 2*N0*fsampling;
            noise = sqrt(NoisePower/2)*(randn(length(signal_hrrc_tx_coded),1)+1i*randn(length(signal_hrrc_tx_coded),1));

            signal_rx = signal_hrrc_tx_coded + noise;
            signal_hhrc_rx = conv(signal_rx, HRC_time);
            signal_hhrc_rx_trunc = signal_hhrc_rx(taps:end-taps+1);

            %% Downsampling
            signal_rx_down = downsample(signal_hhrc_rx_trunc, M);
            
           %% soft coding  
           if decoding==0
            iteration=[0 1 2 3];   %% iteration =0 means without decoding 
            soft_decoded_bits_rx = zeros(Nbits,1);
            for iter=1:4

            for k=1:Npackets
                packet_rx = signal_rx_down(1+(k-1)*codedWordLength : k*codedWordLength); 
                decoded_packet_rx = LdpcSoftDecoder(real((packet_rx)'), H, iteration(iter),N0);
                soft_decoded_bits_rx(1+(k-1)*packetLength:k*packetLength) = decoded_packet_rx(packetLength+1:end);
            end

            BERt(j,iter) = sum(bits_tx ~= soft_decoded_bits_rx)/length(soft_decoded_bits_rx);
            end
           end
        
        
           %% hard decoding
           if decoding==1
            % Demapping
            if(Nbps~=1)
            encoded_bits_rx = (demapping(signal_rx_down,Nbps,modulation))';
            else
            encoded_bits_rx = (demapping(real(signal_rx_down),Nbps,modulation))';   
            end
           iteration=[0 1 2 3];   %% iteration =0 means without decoding 
           for iter=1:4
           decoded_bits_rx = zeros(Nbits,1);
            for k=1:Npackets
                packet_rx = encoded_bits_rx(1+(k-1)*codedWordLength : k*codedWordLength);
                decoded_packet_rx = LdpcHardDecoder(packet_rx, H, iteration(iter));
                decoded_bits_rx(1+(k-1)*packetLength:k*packetLength) = decoded_packet_rx(packetLength+1:end);
            end
           
            BERt(j,iter) = sum(bits_tx ~= decoded_bits_rx)/length(decoded_bits_rx');
          end

           end
        end
BER=BERt;
%% Plot BER results
% figure
SNR = 1:1:20;
plot(SNR,log(BER(:,1)), SNR, log(BER(:,2)),'--', SNR, log(BER(:,3)), '-.', SNR, log(BER(:,4)),'-*')
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('BER curves')
grid on;
legend('PAM','PAM 1 it','PAM 2 it','PAM 3 it');