clc;
clear;
%% Code parameters
% Modulation schemes available using MATLAB's toolbox, which will be used
% in a complex-valued channel
modlist = {'pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM'}; 
bpsList = [1 1 2 4 6 8];
% Pick the modulation
modulation = 'BPSK';
% Determine the number of bits per symbol in the modulation
nmodbits = bpsList(strcmpi(modlist,modulation));
nmax = 256;
kmax = 240;
n = 128;
k = 116;
load('capolar_k_116_n_128_ul.mat')
% n = 127;
% k = 113;
% load('BCH_k_113_n_127.mat')
%% Monte-Carlo parameters
% Sim range in Eb/N0
ebn0 = 4:0.5:7;
% Convert Eb/N0 to SNR
snr_db = ebn0+10*log10(k/n)+10*log10(nmodbits);
numSNR = length(snr_db);
% Pad modulation with meaningless bits if necessary.
mod_pad=zeros(1,mod(n,nmodbits));
for ii=1:numSNR
    % Noise variance
    sigma2 = 1/(10^(0.1*snr_db(ii))); 
    % Using MATLAB's channel function
    awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);
    % Encode
    % Uniform random information word
    u = binornd(ones(1,k),0.5);
    % Codeword
    c = mod(u*G,2);
    % Modulate
    modOut = nrSymbolModulate([c mod_pad]',modulation);
    % Add White Gaussian noise
    rSig = awgnchan(modOut);
    % Soft demodulate
    y_soft = nrSymbolDemodulate(rSig,modulation,sigma2);
    y_soft = y_soft(1:n)';
end
%% Convert LLRs (-LLR_max...+LLR_max) to 6-bit {sign, magnitude} binary string

B = 6;            % total bits (1 sign + 5 magnitude)
LLR = y_soft;     % LLRs
LLR_max = max(abs(y_soft));     % maximum LLR

% Clip to range [-LLR_max, +LLR_max]
LLR = max(min(LLR, LLR_max), -LLR_max);

% Compute magnitude(scaled to 0–31),(close LLRs quantize to the same value)
LLR_mag = round(abs(LLR) / LLR_max * (2^(B-1) - 1));  % 0–31 for 5 bits

% Sign bit: 1 for negative, 0 for positive
y_hard = LLR < 0;
Hy = mod(H*y_hard',2);
% Combine into binary string for each LLR
LLR_bin = strings(1, n);
for i = 1:n
    LLR_mag_bin = dec2bin(LLR_mag(i), B-1);  % 5 bits magnitude
    y_hard_char = num2str(y_hard(i)); % 1 bit sign
    LLR_bin(i) = strcat(y_hard_char, LLR_mag_bin);
end

% 5. Concatenate into one long binary string for VHDL simulation
LLR_concat = strjoin(flip(LLR_bin), '');

% 6. Print result 
fprintf('LLR_in = "%s"\n', LLR_concat);
fprintf("Total bits: %d (expected %d)\n", strlength(LLR_concat), n*B);
%% Streaming LLR input txt for System TCU
LLR_dec = LLR_dec(:)';   % ensure row vector

% Pad with zeros if needed
if length(LLR_dec) < nmax
    LLR_dec = [LLR_dec, zeros(1, nmax - length(LLR_dec))];
end

% Trim if longer (just in case)
LLR_dec = LLR_dec(1:nmax);

% Save to file
dlmwrite('llrs_dec.txt', LLR_dec(:), 'delimiter', '\n');

%% Compare LLR_mag MATLAB and VHDL results
mag_bin_matlab = strings(1, length(LLR_mag));
for i = 1:length(LLR_mag)
    mag_bin_matlab(i) = dec2bin(LLR_mag(i), B-1); % 5 bits each
end
LLR_mag_matlab  = strjoin(flip(mag_bin_matlab), '');

total_bits = nmax * (B - 1);
used_bits  = n * (B - 1);
pad_bits   = total_bits - used_bits;

% Create padding (all zeros)
padding = repmat('0', 1, pad_bits);

% Prepend padding to match VHDL: (nmax*B - 1 downto n*B) <= 0
LLR_mag_matlab_pad = strcat(padding, LLR_mag_matlab);

fprintf("Final bit width: %d (expected %d)\n", strlength(LLR_mag_matlab_pad), total_bits);
fprintf('LLR_mag = "%s"\n', LLR_mag_matlab_pad);
% %%
% LLR_mag_vhdl="00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010111100010101110010011100111001101011110111110101011100011001000100001011010000010100100001101101111000101100011110100010001010111110001110000111000001010011000010101000101000111101111100101100010100011101001110110010110111110011001000110110100100101001001011101010011110000011100110110000001000101110101100011101010011101101010100111110001101110001111000101110001001100011010000101010111110100101100011001111100001010010010010010100001101100111001101111100110101100010010001001101100100000111101111100011110001101010000111001101101111100110000010110100001110010001010010011010111010111001111110101110000101111000001110100110110110001";
% LLR_mag_comp = [LLR_mag_matlab_pad; LLR_mag_vhdl];
% if LLR_mag_comp(1) == LLR_mag_comp(2)
%     disp("The bitstrings are identical!");
% else
%     disp("The bitstrings differ!");
% end
%% Compare y_hard MATLAB and VHDL results
yhard_matlab = strings(1, length(y_hard));
for i = 1:length(y_hard)
    yhard_matlab(i) = num2str(y_hard(i));
end
y_hard_matlab  = strjoin(flip(yhard_matlab), '');

total_bits_yhard = nmax;
used_bits_yhard  = n;
pad_bits_yhard   = total_bits_yhard - used_bits_yhard;

% Create padding (all zeros)
padding_yhard = repmat('0', 1, pad_bits_yhard);

% Prepend padding to match VHDL:
y_hard_matlab = strcat(padding_yhard, y_hard_matlab);

fprintf("Final bit width: %d (expected %d)\n", strlength(y_hard_matlab), total_bits_yhard);
fprintf('y_hard = "%s"\n', y_hard_matlab);
% %%
% y_hard_vhdl = "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011001011011010111100110011001011101000111111111101111101000000011001011011110100000000111011011010111010010000111011011011100";
% y_hard_comp = [y_hard_matlab; y_hard_vhdl];
% if y_hard_comp(1) == y_hard_comp(2)
%     disp("The bitstrings are identical!");
% else
%     disp("The bitstrings differ!");
% end