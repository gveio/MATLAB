% Simulation peforms decoding with OBRGRAND
clear; 
DECODER='ORBGRAND';
%% Code parameters
% Pick the code
% RLC, PAC, CAPOLAR, BCH, eBCH or CRC. 
code_class = 'CRC';

% Random Linear Code
if isequal(code_class,'RLC')
    n=128;
    k=116;
    
% Make a random parity check matrix
[G,H] = make_RLC(k,n,0.5);
code.G=G;
code.H=H;
   
% Polar-assisted convolutional codes
elseif isequal(code_class,'PAC')
    n=128;
    k=116;
    % Generator for convolutional code, taken from arXiv:1908.09594.
    conv_code = [1 1 1 0 1 1 0 1];
    [G, H] = make_pac_code(n,k,conv_code);
    code.G = G;
    code.H = H;
    % Record what convolution was used.
    code.conv_code = conv_code;

% CRC-Assisted Polar code
elseif isequal(code_class,'CAPOLAR')
    % 5G New Radio uplink code
    filename = 'capolar_k_116_n_128_ul.mat';
    code = open(filename);
    G = code.G;
    H = code.H;
    % k,n is the code length
    [k,n]=size(G);

% BCH
elseif isequal(code_class,'BCH')
n           = 127;
k           = 113;
[G, H]      = getGH_BCH(n, k); 
code.G = G;
code.H = H;

% eBCH
%extra parity bit ensures even-weight codewords, reducing error probability
elseif isequal(code_class,'eBCH')
n           = 128;
k           = 113;
[G, H]      = getGH_sys_eBCH(n, k); 
code.G = G;
code.H = H;

% CRC
elseif isequal(code_class,'CRC')
    k=116;
    % Polynomial from https://users.ece.cmu.edu/~koopman/crc/ which has a
    % repository of high quality CRCs curated and maintained by 
    % Philip Koopman, Carnegie Mellon University.
    hex_poly = '0x8f3'; % 12 bit
    poly=koopman2matlab(hex_poly); % Convert from Koopman notation
    % Record the polynomial
    code.poly = poly;
    % Set up the CRC check using MATLAB toolbox
    [G,H,n] = make_CRC_GH(k,poly);
    code.G=G;
    code.H=H;
end

% Record the code_class in the data construct
code.class = code_class;
%% Decoder parameters
Tmax = Inf; % Can reduce to an abandonment value.
%% Monte-Carlo parameters
% Sim range in Eb/N0
EbN0dB = 4:0.5:6;
% Run simulation until this many errors observed for each Eb/N0 value
NoErrors = 50;
maxIt       = 10^6; % Bounds on simulation iterations
minIt       = 10^2;
%% Code and channel
R           = (k / n);
EsN0dB      = EbN0dB + 10 * log10(2*R);
numSNR      = length(EsN0dB);
%% Loop over SNRs
BLER        = zeros(1, numSNR);
BER         = zeros(1, numSNR);
NGavg       = zeros(1, numSNR);
NGavg_p     = zeros(1, numSNR);
for ii = 1:numSNR
    BlockError  = 0;
    BitError    = 0;
    ntx         = 0;
    NG          = 0;
    sigma = 1 / sqrt(10^(EsN0dB(ii) / 10));
    while ((BlockError < NoErrors && ntx < maxIt) || ntx < minIt)
        ntx = ntx + 1;
        c = zeros(1, n); % To store codeword.
        %% Encode
        u = randsrc(1, k, [0, 1]); % Random information bits
        c= mod(u*G, 2);
        %% binary input AWGN channel
        x = 1 - 2 * c; % BPSK mapping
        y = x + sigma * randn([1, n]); % AWGN channel
        L_channel = 2 * y / (sigma^2); % LLR
        %% Decoding
        [y_decoded, n_guess] = ORBGRAND_mex(L_channel',uint8(reshape(H', [], 1)),uint64(Tmax));
        y_decoded =y_decoded';
        NG = NG + n_guess;
        %% error collection
        if (~isequal(c, y_decoded))  % If decoded word differs from transmitted
            BlockError = BlockError + 1; % Count as a block error
            %BitError = BitError + sum(abs(y_decoded-c));
            BitError = BitError + sum(y_decoded(:) ~= c(:));
            % uhat = y_decoded(1, 1:k);
            % BitError = BitError + sum(uhat(:) ~= u(:));
            % use this metric when comparing the 
            % decoded message bits (uhat) to the original message bits (u).
        end     
    end
    disp(['---' code_class '[' num2str(n) ',' num2str(k) ']---Eb/N0 dB ', num2str(EbN0dB(ii)), ' dB:---'])
    BLER    = BlockError / ntx;
    BER     = BitError / (ntx * n);
    % BER(ii)     = BitError / (ntx * k); 
    % use this metric when comparing the decoded message bits (uhat) 
    % to the original message bits (u).
    NGavg   = NG / (ntx);
    disp([' Blocks decoded   = ', num2str(ntx)]);
    disp([' BLER             = ', num2str(BLER)]);
    disp([' BER              = ', num2str(BER)]);
    disp([' NGavg            = ', num2str(NGavg)]);    
    % disp([' NGavg/(info bit) = ', num2str(NGavg(ii)/k)]);
end
    % save(['../RESULTS_sim/ORBGRAND -' code_class '-' num2str(n) '-' num2str(k) '.mat'], 'EbN0dB', 'BLER', 'BER', 'NGavg', 'NGavg_p', 'G', 'H');
