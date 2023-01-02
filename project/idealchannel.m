% %% Projet modulation & coding
addpath(genpath('function'));
%% Parameters
clear all
f_cut = 1e+6; % cut off frequency of the nyquist filter [Mhz]
M = 10; % oversampling factor
fsymb = 2*f_cut; % the symbol frequency 
fsampling = M*fsymb; % sampling frequency
Tsymb = 1/fsymb; % time between two symbols.
Npbsa=[1 2 4 6];
beta = 0.3; % roll-off factor
Npackets = 6; % choose such that Nbits/Nbps is an integer
packetLength = 128;
codedWordLength = 2*packetLength; % make sure rate=1/2  bits_tex/codeword = 1/2
Nbits = Npackets*packetLength; % bit stream length
NcodedBits = Npackets*codedWordLength; % full coded word length

bits_tx = randi(2,Nbits,1)-1; % generate random 0101110011,,,,,, 
bits_tx_coded = zeros(NcodedBits,1);
%% main loop
for jj=1:4
    

    Nbps=Npbsa(jj);
    %% chanl coding      
    H0 = makeLdpc(packetLength, codedWordLength, 0, 1, 3);
for k=1:Npackets
    packet_tx = bits_tx(1+(k-1)*packetLength : k*packetLength);% choose each packet bits_tex value , each packet length=128
    [codedbits, H] = makeParityChk(packet_tx , H0, 0);
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
    EB_N0 = 1:0.2:20;
    BER1 = zeros(length(EB_N0),5);
    symnol_energy_coded = (trapz(abs(signal_hrrc_tx_coded).^2))/fsampling; % total symbol energy
    Eb = symnol_energy_coded/Nbits; % total symbol energy / total bits = bit energy
    Eb =Eb/2;%%% pass band bit energy = baseband bit energy/2

%Bit rate = Number of bits*Sampling frequency s  
%Eb Calculated as received signal power (in Watts) divided by bit rate (in 1s):
        for j = 1:length(EB_N0)
            N0 = Eb/(10.^(EB_N0(j)/10));%%noise energy= ((N0/2)w*2)*T=N0  //w=1/T T is symbol period 
            Noisepower= N0*fsampling; 
            % psd of baseband=N0 ; psd of passband=N0/2
            % bandwidth of baseband =fsampling   passband=  positive + negative = fsampling +fsampling
            % Noisepower of baseband : noisepower_base_real=N0*fsampling/2   noisepower_base_imaginery=N0*fsampling/2   
            % Pass band_noisepower_(psd of passband=N0/2):  positive freuqency ：(N0/2)*fsampling ; negative freuqency： (N0/2)*fsampling 
            % due to wss noise, the mean value of wss noise is 0 => the value of noisepower is equal to it's variance  
            %%https://dsp.stackexchange.com/questions/46385/variance-in-the-time-domain-versus-variance-in-frequency-domain
            Noise_var=Noisepower;
            noise = sqrt(Noise_var)*(randn(length(signal_hrrc_tx_coded),1)+1i*randn(length(signal_hrrc_tx_coded),1));
            % Noisepower = noise variance =>  std= sqrt(Noisepower)
            signal_rx = signal_hrrc_tx_coded + noise;
            signal_hhrc_rx = conv(signal_rx, HRC_time);
            signal_hhrc_rx_trunc = signal_hhrc_rx(taps:end-taps+1);

            %% Downsampling
            signal_rx_down = downsample(signal_hhrc_rx_trunc, M);
                   

            %% Demapping
            if(Nbps~=1)
            encoded_bits_rx = (demapping(signal_rx_down,Nbps,modulation))';
            else
            encoded_bits_rx = (demapping(real(signal_rx_down),Nbps,modulation))';   
            end
           %% hard decoding
           decoded_bits_rx = zeros(Nbits,1);
            for k=1:Npackets
                packet_rx = encoded_bits_rx(1+(k-1)*codedWordLength : k*codedWordLength);
                decoded_packet_rx = LdpcHardDecoder(packet_rx, H, 5);
                decoded_bits_rx(1+(k-1)*packetLength:k*packetLength) = decoded_packet_rx(packetLength+1:end);
            end
           
            BER(j,jj) = sum(bits_tx ~= decoded_bits_rx)/length(decoded_bits_rx');
        end
      
        
end
%% plot psd
figure
pwelch(signal_rx)
hold on
pwelch(signal_hhrc_rx_trunc)
%% plot HRC AND RC impulse response 
t=-182*Tsymb:Tsymb:182*Tsymb;
figure
plot(t,normal_Rc_time);
hold on
plot(t,HRC_time);
legend('Raised cosin filter','half-root Raised cosin filter')
title('Riased cosin filter in time domain')
xlabel('time(s)')
ylabel('Normalized Amplitude')
%% plot HRC AND RC in freq domain
figure 
plot(f,normal_Rc_freq)
hold on
plot(f,HRC_freq)
xlabel('f')
ylabel('Amplitude')
title('Riased cosin filter in freq domain')


%% Plot BER results
figure
for jj=1:4
plot(EB_N0,log(BER(:,jj)))
xlabel('E_B/N_0 [dB]');
ylabel('BER');
legend('PAM','4QAM','16QAM','64QAM');
title('BER curves')
grid on;
hold on
end