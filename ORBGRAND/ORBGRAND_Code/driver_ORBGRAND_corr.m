% ==========================================================
% ORBGRAND FER sweep under AR(1) correlated noise
%
% Includes:
%   1) Classic AWGN LLR
%   2) Correlation-aware LLR from Eq. (15)
%   3) Sign-convention sanity check
%   4) LLR magnitude sanity check
%
% Channel:
%   z(i) = rho*z(i-1) + sqrt(1-rho^2)*w(i)
%   rSig = modOut + sqrt(sigma2)*z
%
% rho = 0 gives standard AWGN
% ==========================================================

clear; clc;

%% Decoder selection
DECODER = 'ORBGRAND-BASELINE';
% DECODER = 'ORBGRAND-MSB-3';
% DECODER = 'ORBGRAND_C';
% DECODER = 'ORBGRAND_M';

%% LLR mode
% 'classic'  : classic AWGN LLR
% 'corr_eq15': correlation-aware LLR from Eq. (15)
LLR_MODE = 'corr_eq15';

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
ebn0 = 4:0.5:9;
err_thresh = 50;

rho_list = [0 0.5 0.8];

snr_db = ebn0 + 10*log10(R) + 10*log10(nmodbits);

numSNR = length(snr_db);
numRho = length(rho_list);

mod_pad = zeros(1,mod(n,nmodbits));

%% Accumulators: rows = rho, cols = Eb/N0
num_decoded        = zeros(numRho,numSNR);
num_demod_errs     = zeros(numRho,numSNR);
num_demod_bit_errs = zeros(numRho,numSNR);
num_errs           = zeros(numRho,numSNR);
num_bit_errs       = zeros(numRho,numSNR);
num_queries        = zeros(numRho,numSNR);

%% Extra sanity statistics
mean_abs_llr_classic = zeros(numRho,numSNR);
mean_abs_llr_corr    = zeros(numRho,numSNR);
hd_ber_classic       = zeros(numRho,numSNR);
hd_ber_corr_flip     = zeros(numRho,numSNR);
hd_ber_corr_noflip   = zeros(numRho,numSNR);

count_sanity         = zeros(numRho,numSNR);

%% Main simulation
for rr = 1:numRho

    rho = rho_list(rr);

    fprintf('\n========================================\n');
    fprintf('Running rho = %.2f\n', rho);
    fprintf('LLR mode = %s\n', LLR_MODE);
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
            modOut = real(modOut(1:n).');

            %% AR(1)-correlated Gaussian noise
            w = randn(1,n);
            z = zeros(1,n);
            z(1) = w(1);

            for jj = 2:n
                z(jj) = rho*z(jj-1) + sqrt(1-rho^2)*w(jj);
            end

            noise = sqrt(sigma2) * z;

            %% Received signal
            rSig = modOut + noise;
            y = real(rSig(:)).';   % row vector

            %% Classic LLR
            % MATLAB/ORBGRAND convention used here:
            % positive LLR -> bit 0 / BPSK +1 is more likely
            y_soft_classic = 2*y/sigma2;

            %% Correlation-aware LLR from Eq. (15)
            y_prev = [0, y(1:end-1)];

            if abs(rho) < 1e-12
                lambda_paper = -2*y/sigma2;
            else
                lambda_paper = ( ...
                    -2*y + 2*rho*y_prev ...
                    + abs(y_prev - rho*y - rho) ...
                    - abs(y_prev - rho*y + rho) ...
                    ) ./ (sigma2*(1-rho^2));
            end

            % Two possible sign conventions:
            y_soft_corr_noflip = lambda_paper;
            y_soft_corr_flip   = -lambda_paper;

            %% Sanity checks: hard decisions
            hd_classic = (y_soft_classic < 0);
            hd_corr_flip = (y_soft_corr_flip < 0);
            hd_corr_noflip = (y_soft_corr_noflip < 0);

            count_sanity(rr,ii) = count_sanity(rr,ii) + 1;

            mean_abs_llr_classic(rr,ii) = mean_abs_llr_classic(rr,ii) + mean(abs(y_soft_classic));
            mean_abs_llr_corr(rr,ii)    = mean_abs_llr_corr(rr,ii)    + mean(abs(y_soft_corr_flip));

            hd_ber_classic(rr,ii)     = hd_ber_classic(rr,ii)     + mean(hd_classic ~= c);
            hd_ber_corr_flip(rr,ii)   = hd_ber_corr_flip(rr,ii)   + mean(hd_corr_flip ~= c);
            hd_ber_corr_noflip(rr,ii) = hd_ber_corr_noflip(rr,ii) + mean(hd_corr_noflip ~= c);

            %% Select LLR for decoding
            switch lower(LLR_MODE)

                case 'classic'
                    y_soft = y_soft_classic;

                case 'corr_eq15'
                    % Use the sign-flipped version because this matches
                    % the classic convention at rho=0.
                    y_soft = y_soft_corr_flip;

                otherwise
                    error('Unsupported LLR_MODE.');
            end

            y_soft(1) = 2*y(1)/sigma2;

            y_soft = double(real(y_soft(:))).';   % force real double row

            %% Hard demodulate using selected LLR
            y_demod = (y_soft < 0);

            if ~isequal(y_demod, c)
                num_demod_errs(rr,ii) = num_demod_errs(rr,ii) + 1;
                num_demod_bit_errs(rr,ii) = num_demod_bit_errs(rr,ii) + sum(abs(y_demod-c));
            end

            %% Decode
            if isequal(DECODER,'ORBGRAND-MSB-3')

                [y_decoded,n_guess] = ORBGRAND_msb3_quant_mex( ...
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

        %% Print sanity summary for this point
        mean_abs_llr_classic(rr,ii) = mean_abs_llr_classic(rr,ii) / count_sanity(rr,ii);
        mean_abs_llr_corr(rr,ii)    = mean_abs_llr_corr(rr,ii)    / count_sanity(rr,ii);

        hd_ber_classic(rr,ii)     = hd_ber_classic(rr,ii)     / count_sanity(rr,ii);
        hd_ber_corr_flip(rr,ii)   = hd_ber_corr_flip(rr,ii)   / count_sanity(rr,ii);
        hd_ber_corr_noflip(rr,ii) = hd_ber_corr_noflip(rr,ii) / count_sanity(rr,ii);

        fprintf('\nSANITY CHECK rho=%.2f, Eb/N0=%.2f dB:\n', rho, ebn0(ii));
        fprintf('  mean|LLR classic|       = %.4f\n', mean_abs_llr_classic(rr,ii));
        fprintf('  mean|LLR corr Eq15|     = %.4f\n', mean_abs_llr_corr(rr,ii));
        fprintf('  HD BER classic          = %.4e\n', hd_ber_classic(rr,ii));
        fprintf('  HD BER corr flip        = %.4e\n', hd_ber_corr_flip(rr,ii));
        fprintf('  HD BER corr no-flip     = %.4e\n', hd_ber_corr_noflip(rr,ii));

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
result.LLR_MODE = LLR_MODE;
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

result.mean_abs_llr_classic = mean_abs_llr_classic;
result.mean_abs_llr_corr = mean_abs_llr_corr;
result.hd_ber_classic = hd_ber_classic;
result.hd_ber_corr_flip = hd_ber_corr_flip;
result.hd_ber_corr_noflip = hd_ber_corr_noflip;

result.G = G;
result.H = H;

if isequal(code_class,'CRC')
    filename = ['../RESULTS/AR1_FER_SWEEP_' LLR_MODE '_' DECODER '_' code_class '_' hex_poly '_' ...
        num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
else
    filename = ['../RESULTS/AR1_FER_SWEEP_' LLR_MODE '_' DECODER '_' code_class '_' ...
        num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
end

save(filename,'result');

fprintf('\nSaved results to:\n%s\n', filename);

%% FER plot
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
title(sprintf('%s, %s (%d,%d), %s', ...
    DECODER, code_class, n, k, LLR_MODE), 'Interpreter','none');

%% BER plot
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
title(sprintf('%s, %s (%d,%d), %s', ...
    DECODER, code_class, n, k, LLR_MODE), 'Interpreter','none');

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
title(sprintf('%s, %s (%d,%d), query complexity, %s', ...
    DECODER, code_class, n, k, LLR_MODE), 'Interpreter','none');

%% LLR magnitude sanity plot
figure; clf;
for rr = 1:numRho
    plot(ebn0, mean_abs_llr_corr(rr,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('corr LLR, \\rho = %.1f', rho_list(rr)));
    hold on;
end
plot(ebn0, mean_abs_llr_classic(1,:), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'classic LLR reference');
grid on;
xlabel('$E_b/N_0$ (dB)','Interpreter','latex');
ylabel('Mean $|\mathrm{LLR}|$','Interpreter','latex');
legend('Location','northwest');
title('LLR magnitude sanity check','Interpreter','latex');

%% Hard-decision BER sanity plot
figure; clf;
for rr = 1:numRho
    semilogy(ebn0, hd_ber_corr_flip(rr,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('corr flip, \\rho = %.1f', rho_list(rr)));
    hold on;
end
semilogy(ebn0, hd_ber_classic(1,:), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'classic reference');
grid on;
xlabel('$E_b/N_0$ (dB)','Interpreter','latex');
ylabel('Hard-decision BER','Interpreter','latex');
legend('Location','southwest');
title('Hard-decision BER sanity check','Interpreter','latex');