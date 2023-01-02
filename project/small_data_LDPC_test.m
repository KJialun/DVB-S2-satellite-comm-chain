
%% init
Npackets = 1; % choose such that Nbits/Nbps is an integer
packetLength = 128;
codedWordLength = 2*packetLength; % make sure rate=1/2  bits_tex/codeword = 1/2
Nbits = Npackets*packetLength; % bit stream length
NcodedBits = Npackets*codedWordLength; % full coded word length

bits_tx = randi(2,Nbits,1)-1; % generate random 0101110011,,,,,, 
bits_tx_coded = zeros(NcodedBits,1);

%% encoding 
H0 = makeLdpc(packetLength, codedWordLength, 0, 1, 3);
for k=1:Npackets
    packet_tx = bits_tx(1+(k-1)*packetLength : k*packetLength);% choose each packet bits_tex value , each packet length=128
    [codedbits, H] = makeParityChk(packet_tx , H0, 0);
    bits_tx_coded(1+(k-1)*codedWordLength : k*codedWordLength) = [codedbits packet_tx];
end
%% MAPPING CHANNEL ..... 

encoded_bits_rx=bits_tx_coded';

%% Hard Decoding
decoded_bits_rx = zeros(Nbits,1);
for k=1:Npackets
    packet_rx = encoded_bits_rx(1+(k-1)*codedWordLength : k*codedWordLength);
    decoded_packet_rx = LdpcHardDecoder(packet_rx, H, 2);
    decoded_bits_rx(1+(k-1)*packetLength:k*packetLength) = decoded_packet_rx(packetLength+1:end);
end






