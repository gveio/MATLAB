% ==========================================================
% ORBGRAND FER sweep under AR(1) correlated noise
%
% Channel model:
%   z(n) = rho*z(n-1) + sqrt(1-rho^2)*w(n)
%   rSig = modOut + sqrt(sigma2)*z
%
% rho = 0 gives standard AWGN
% ==========================================================

clear; clc;

% DECODER = 'ORBGRAND-MSB3-TIE-POLICY';
DECODER = 'ORBGRAND-BASELINE';
% DECODER = 'ORBGRAND_C';
% DECODER = 'ORBGRAND_M';

%% Code parameters
modlist = {'pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM'}; 
bpsList = [1 1 2 4 6 8];

modulation = 'BPSK';
nmodbits = bpsList(strcmpi(modlist,modulation));

code_class = 'CRC';

switch upper(code_class)

    case 'CAPOLAR'
        filename = 'capolar_k_116_n_128_ul.mat';
        code = open(filename);
        G = code.G;
        H = code.H;
        [k,n] = size(G);

    case 'CRC'
        k = 240;
        hex_poly = '0xd175';
        poly = koopman2matlab(hex_poly);
        [G,H,n] = make_CRC_GH(k,poly);

    otherwise
        error('Unsupported code_class');
end

R = k/n;

%% Decoder parameters
max_query = Inf;

%% Monte Carlo parameters
ebn0 = 4:0.5:7;
err_thresh = 20;

rho_list = [0 0.5 0.8];

snr_db = ebn0 + 10*log10(R) + 10*log10(nmodbits);

numSNR = length(snr_db);
numRho = length(rho_list);

mod_pad = zeros(1,mod(n,nmodbits));

%% Accumulators: rows = rho, columns = Eb/N0
num_decoded        = zeros(numRho,numSNR);
num_demod_errs     = zeros(numRho,numSNR);
num_demod_bit_errs = zeros(numRho,numSNR);
num_errs           = zeros(numRho,numSNR);
num_bit_errs       = zeros(numRho,numSNR);
num_queries        = zeros(numRho,numSNR);

%% Main simulation
for rr = 1:numRho

    rho = rho_list(rr);

    fprintf('\n========================================\n');
    fprintf('Running rho = %.2f\n', rho);
    fprintf('========================================\n');

    for ii = 1:numSNR

        sigma2 = 1/(10^(0.1*snr_db(ii)));

        fprintf('\nEb/N0 = %.2f dB, SNR = %.4f dB, sigma2 = %.6g\n', ...
            ebn0(ii), snr_db(ii), sigma2);

        while num_errs(rr,ii) < err_thresh

            num_decoded(rr,ii) = num_decoded(rr,ii) + 1;

            %% Encode
            u = binornd(ones(1,k),0.5);
            c = mod(u*G,2);

            %% Modulate
            modOut = nrSymbolModulate([c mod_pad]', modulation);
            modOut = modOut(1:n).';

            %% AR(1) correlated Gaussian noise
            % z has unit variance.
            % rho = 0 gives standard independent AWGN.
            w = randn(1,n);
            z = zeros(1,n);
            z(1) = w(1);

            for jj = 2:n
                z(jj) = rho*z(jj-1) + sqrt(1-rho^2)*w(jj);
            end

            noise = sqrt(sigma2) * z;

            %% Received signal
            rSig = modOut + noise;

            %% Soft demodulate
            y_soft = nrSymbolDemodulate(rSig.', modulation, sigma2);
            y_soft = y_soft(1:n).';

            %% Hard demodulate
            y_demod = (y_soft < 0);

            if ~isequal(y_demod, c)
                num_demod_errs(rr,ii) = num_demod_errs(rr,ii) + 1;
                num_demod_bit_errs(rr,ii) = num_demod_bit_errs(rr,ii) + sum(abs(y_demod-c));
            end

            %% Decode
            if isequal(DECODER,'ORBGRAND-MSB3-TIE-POLICY')

                [y_decoded,n_guess] = ORBGRAND_msb3_tie_policy_quant_mex( ...
                    y_soft', uint8(reshape(H', [], 1)), uint64(max_query));
                y_decoded = y_decoded';

            elseif isequal(DECODER,'ORBGRAND_C')

                [y_decoded,n_guess] = ORBGRAND_mex( ...
                    y_soft', uint8(reshape(H', [], 1)), uint64(max_query));
                y_decoded = y_decoded';

            elseif isequal(DECODER,'ORBGRAND-BASELINE')

                [y_decoded,n_guess] = ORBGRAND_quant_mex( ...
                    y_soft', uint8(reshape(H', [], 1)), uint64(max_query));
                y_decoded = y_decoded';

            elseif isequal(DECODER,'ORBGRAND_M')

                [y_decoded,~,n_guess,~] = bin_ORBGRAND(H,max_query,y_soft); 

            else
                error('DECODER not found');
            end

            num_queries(rr,ii) = num_queries(rr,ii) + n_guess;

            %% Decode error
            if ~isequal(y_decoded, c)

                num_errs(rr,ii) = num_errs(rr,ii) + 1;
                num_bit_errs(rr,ii) = num_bit_errs(rr,ii) + sum(abs(y_decoded-c));

                fprintf('rho=%.2f, Eb/N0=%.2f dB, Dec=%d, Err=%d, FER=%.4e, BER=%.4e, EG=%.0f\n', ...
                    rho, ...
                    ebn0(ii), ...
                    num_decoded(rr,ii), ...
                    num_errs(rr,ii), ...
                    num_errs(rr,ii)/num_decoded(rr,ii), ...
                    num_bit_errs(rr,ii)/(num_decoded(rr,ii)*n), ...
                    num_queries(rr,ii)/num_decoded(rr,ii));
            end
        end
    end
end

%% Summary statistics
BLERdemod = num_demod_errs ./ num_decoded;
BERdemod  = num_demod_bit_errs ./ (num_decoded*n);

FER = num_errs ./ num_decoded;
BER = num_bit_errs ./ (num_decoded*n);
EG  = num_queries ./ num_decoded;

%% Store results
result = struct();

result.DECODER = DECODER;
result.code_class = code_class;
result.n = n;
result.k = k;
result.R = R;
result.modulation = modulation;
result.nmodbits = nmodbits;

result.ebn0 = ebn0;
result.snr_db = snr_db;
result.rho_list = rho_list;
result.max_query = max_query;
result.err_thresh = err_thresh;

result.num_decoded = num_decoded;
result.num_demod_errs = num_demod_errs;
result.num_demod_bit_errs = num_demod_bit_errs;
result.num_errs = num_errs;
result.num_bit_errs = num_bit_errs;
result.num_queries = num_queries;

result.BLERdemod = BLERdemod;
result.BERdemod = BERdemod;
result.FER = FER;
result.BER = BER;
result.EG = EG;

result.G = G;
result.H = H;

if isequal(code_class,'CRC')
    filename = ['../RESULTS/AR1_FER_SWEEP_' DECODER '_' code_class '_' hex_poly '_' ...
        num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
else
    filename = ['../RESULTS/AR1_FER_SWEEP_' DECODER '_' code_class '_' ...
        num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
end

save(filename,'result');

fprintf('\nSaved results to:\n%s\n', filename);

%% Quick FER plot
figure; clf;
for rr = 1:numRho
    semilogy(ebn0, FER(rr,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\rho = %.1f', rho_list(rr)));
    hold on;
end
grid on;
xlabel('$E_b/N_0$ (dB)','Interpreter','latex');
ylabel('FER','Interpreter','latex');
legend('Location','southwest');
title(sprintf('%s, %s (%d,%d), AR(1) correlated noise', ...
    DECODER, code_class, n, k), 'Interpreter','none');

%% Quick BER plot
figure; clf;
for rr = 1:numRho
    semilogy(ebn0, BER(rr,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\rho = %.1f', rho_list(rr)));
    hold on;
end
grid on;
xlabel('$E_b/N_0$ (dB)','Interpreter','latex');
ylabel('BER','Interpreter','latex');
legend('Location','southwest');
title(sprintf('%s, %s (%d,%d), AR(1) correlated noise', ...
    DECODER, code_class, n, k), 'Interpreter','none');

%% Query complexity plot
figure; clf;
for rr = 1:numRho
    semilogy(ebn0, EG(rr,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\rho = %.1f', rho_list(rr)));
    hold on;
end
grid on;
xlabel('$E_b/N_0$ (dB)','Interpreter','latex');
ylabel('Avg. codebook queries','Interpreter','latex');
legend('Location','northeast');
title(sprintf('%s, %s (%d,%d), query complexity', ...
    DECODER, code_class, n, k), 'Interpreter','none');