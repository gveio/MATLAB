% Simulation peforms decoding with OBRGRAND
clear;
DECODER='ORBGRAND-MSB-4';
%% Code parameters
% Modulation schemes available using MATLAB's toolbox, which will be used
% in a complex-valued channel
modlist = {'pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM'}; 
bpsList = [1 1 2 4 6 8];
% Pick the modulation
modulation = 'BPSK';
% Determine the number of bits per symbol in the modulation
nmodbits = bpsList(strcmpi(modlist,modulation));

% Pick the code
% RLC, PAC, CAPOLAR, BCH, eBCH or CRC. 
code_class = 'CAPOLAR';

% Random Linear Code
if isequal(code_class,'RLC')
    n=256;
    k=240;
    
% Make a random parity check matrix
[G,H] = make_RLC(k,n,0.5);
code.G=G;
code.H=H;
   
% Polar-assisted convolutional codes
elseif isequal(code_class,'PAC')
    n=256;
    k=240;
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
    
    % filename = 'capolar_k_240_n_256_ul.mat';
    % code = load(filename);          % load instead of open (recommended)
    % S = code.code;            % struct name
    % G = S.G;
    % H = S.H;
    % [k,n] = size(G);

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
    % k=116;
    k=240;
    % Polynomial from https://users.ece.cmu.edu/~koopman/crc/ which has a
    % repository of high quality CRCs curated and maintained by 
    % Philip Koopman, Carnegie Mellon University.
    % hex_poly = '0x8f3'; % 12 bit
    hex_poly = '0xd175';  % CRC-16
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
max_query = Inf; % Can reduce to an abandonment value.
%% Monte-Carlo parameters
% Sim range in Eb/N0
ebn0 = 4:0.5:7;
% Run simulation until this many errors observed for each Eb/N0 value
err_thresh = 50;
%% Simulation
% Convert Eb/N0 to SNR
snr_db = ebn0+10*log10(k/n)+10*log10(nmodbits);
numSNR      = length(snr_db);
% Pad modulation with meaningless bits if necessary.
mod_pad=zeros(1,mod(n,nmodbits));

%For each SNR record the following
num_decoded = zeros(1,length(snr_db)); % Number of packets decoded
num_demod_errs = num_decoded; % Number of demodulations that are in error
num_demod_bit_errs = num_decoded; % Number of erroneous demodulated bits
num_errs = num_decoded; % Number of erroneous decodings
num_bit_errs = num_decoded; % Number of erroneous bits
num_queries = num_decoded; % Total number of code-book queries

for ii=1:numSNR
    % Noise variance
    sigma2 = 1/(10^(0.1*snr_db(ii))); 
    % Using MATLAB's channel function
    awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);

    % Keep decoding pacekts until you see err_thresh erroneous decodings
    while num_errs(ii)<err_thresh
        num_decoded(ii)=num_decoded(ii)+1;
                
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
        % Hard demodulate
        y_demod = (y_soft<0);

        % Count the times the demodulation is in error
        if ~isequal(y_demod, c)
            num_demod_errs(ii) = num_demod_errs(ii)+1;
            % Count the demodulated bit errors
            num_demod_bit_errs(ii) = num_demod_bit_errs(ii)+sum(abs(y_demod-c));
        end
        if isequal(DECODER,'ORBGRAND-MSB-4')
            [y_decoded,n_guess] = ORBGRAND_msb4_quant_mex(y_soft',uint8(reshape(H', [], 1)),uint64(max_query));
            y_decoded =y_decoded';
        elseif isequal(DECODER,'ORBGRAND_C')
        % Decode with basic ORBGRAND C implementation
            [y_decoded,n_guess] = ORBGRAND_mex(y_soft',uint8(reshape(H', [], 1)),uint64(max_query));
            y_decoded =y_decoded';
        elseif isequal(DECODER,'ORBGRAND-BASELINE')
            % Decode with basic ORBGRAND C implementation with quantization
            [y_decoded,n_guess] = ORBGRAND_quant_mex(y_soft',uint8(reshape(H', [], 1)),uint64(max_query));
            y_decoded =y_decoded';
        elseif isequal(DECODER,'ORBGRAND_M')
            % Decode with basic ORBGRAND MATLAB implementation
            [y_decoded,~,n_guess,~] = bin_ORBGRAND(H,max_query,y_soft); 
        else
            disp('DECODER not found');
            return;
        end
       
        % Total number of queries made at this SNR
        num_queries(ii) = num_queries(ii) + n_guess;
        
        % If there is an error in the decoding
         % If there is an error in the decoding
        if ~isequal(y_decoded, c)
            % Increment the number of errors observed at this SNR
            num_errs(ii)=num_errs(ii)+1;
            % Count the bit errors.
            num_bit_errs(ii) = num_bit_errs(ii) + sum(abs(y_decoded-c));
            % num_bit_errs(ii) = num_bit_errs(ii) + sum(y_decoded(:) ~= c(:));
            % uhat = y_decoded(1, 1:k);
            % num_bit_errs(ii) = num_bit_errs(ii) + sum(uhat(:) ~= u(:));
            % use this metric when comparing the 
            % decoded message bits (uhat) to the original message bits (u).

         % Report to the terminal every time a new error is observed
         disp(['n=' num2str(n) ', k=' num2str(k) ', R=' num2str(k/n,'%.2f') ', ' num2str(ebn0(ii)) ' dB, Dec=' num2str(num_decoded(ii)) ', Err=' num2str(num_errs(ii)) ', BLER=' num2str(num_errs(ii)/num_decoded(ii)) ', BER=' num2str(num_bit_errs(ii)/(num_decoded(ii)*n)) ', EG=' num2str(num_queries(ii)/num_decoded(ii),'%.0f')]);
        
         %disp(['n=' num2str(n) ', k=' num2str(k) ', R=' num2str(k/n,'%.2f') ', ' num2str(ebn0(ii)) ' dB, Dec=' num2str(num_decoded(ii)) ', Err=' num2str(num_errs(ii)) ', BLER=' num2str(num_errs(ii)/num_decoded(ii)) ', BER=' num2str(num_bit_errs(ii)/(num_decoded(ii)*k)) ', EG=' num2str(num_queries(ii)/num_decoded(ii),'%.0f')]);
         % use this metric when comparing the decoded message bits (uhat) 
         % to the original message bits (u)
        end
    end
end

% Save the results.

if isequal(code.class,'CRC')
    filename = ['../RESULTS/' DECODER '_' code.class '_' hex_poly '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
elseif isequal(code.class,'PAC')
    filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(bin2dec(num2str(conv_code))) '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
else
    filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
end

% Store the data in a way where results from future simulations with the 
% same code can be appended.

if exist(filename, 'file') == 2
    load(filename,'code')
    for ii=1:numSNR
        these = find(code.snr==snr_db(ii));
        if isempty(these)
            code.snr(end+1) = snr_db(ii);
            code.num_decoded(end+1)=num_decoded(ii);
            code.num_demod_errs(end+1)=num_demod_errs(ii);
            code.num_demod_bit_errs(end+1)=num_demod_bit_errs(ii);
            code.num_errs(end+1)=num_errs(ii);
            code.num_bit_errs(end+1)=num_bit_errs(ii);
            code.num_queries(end+1)=num_queries(ii);
        else
            code.num_decoded(these)=code.num_decoded(these)+num_decoded(ii);
            code.num_demod_errs(these)=code.num_demod_errs(these)+num_demod_errs(ii);
            code.num_demod_bit_errs(these)=code.num_demod_bit_errs(these)+num_demod_bit_errs(ii);
            code.num_errs(these)=code.num_errs(these)+num_errs(ii);
            code.num_bit_errs(these)=code.num_bit_errs(these)+num_bit_errs(ii);
            code.num_queries(these)=code.num_queries(these)+num_queries(ii);
        end
    end     
else
    code.n = n;
    code.k = k;
    code.modulation = modulation;
    code.nmodbits = nmodbits;
    code.G = G;
    code.snr = snr_db;
    code.num_decoded=num_decoded;
    code.num_demod_errs=num_demod_errs;
    code.num_demod_bit_errs=num_demod_bit_errs;
    code.num_errs=num_errs;
    code.num_bit_errs=num_bit_errs;
    code.num_queries=num_queries;
end

% Sort according to increasing snr;
[~,ord]=sort(code.snr);
code.snr = code.snr(ord);
code.num_decoded = code.num_decoded(ord);
code.num_demod_errs = code.num_demod_errs(ord);
code.num_demod_bit_errs = code.num_demod_bit_errs(ord);
code.num_errs = code.num_errs(ord);
code.num_bit_errs = code.num_bit_errs(ord);
code.num_queries = code.num_queries(ord);

% Summary statistics   
code.ebn0 = code.snr-10*log10(k/n)-10*log10(nmodbits); %Eb/N0
code.R = code.k/code.n; % Code rate
code.BLERdemod = code.num_demod_errs./code.num_decoded; % Demod BLER
code.BERdemod = code.num_demod_bit_errs./(code.num_decoded*code.n); % Demod BER
code.BLER = code.num_errs./code.num_decoded; % Decoded BLER
code.BER = code.num_bit_errs./(code.num_decoded*code.n); % Decoded BER
code.EG = code.num_queries./code.num_decoded; % Average number of code-book queries
code.max_query = max_query;

% Code
code.G=G;
code.H=H;

save(filename,'code')

disp(['Decoder=' DECODER ', code=' code.class ', ' num2str(sum(num_decoded)) ' packets decoded.'])

