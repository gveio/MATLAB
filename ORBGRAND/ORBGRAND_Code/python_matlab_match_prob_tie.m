% ==========================================================
% Monte Carlo framework for step-wise MSB tie probability
% using the precision-aware bitonic sorter MEX
%
% Supported configurations:
%   1) CA-Polar(128,116)
%   2) CRC(256,240)
%
% MEX required:
%   [sorted_idx, tie_count, total_count] = pa_step_mex(y_soft)
%
% Extra output added here:
%   result.p_msb3_mc : empirical MSB3 input distribution per SNR
%
% Notes:
%   - Tie probability at step t:
%         p(t) = (# MSB ties at step t) / (# total valid comparisons at step t)
%   - One step corresponds to one internal bitonic compare-exchange pass.
%   - We also measure the empirical histogram of the 3-MSB magnitude values
%     at the sorter input, to compare directly against the Python model.
% ==========================================================

clear; clc;

%% User options
code_class = 'CAPOLAR';   % 'CAPOLAR' or 'CRC'
ebn0       = 4:1:7;       % Eb/N0 sweep
nFrames    = 2e4;         % Monte Carlo frames per SNR

%% Modulation: BPSK only
modulation = 'BPSK';
nmodbits   = 1;

%% Load selected code
switch upper(code_class)

    case 'CAPOLAR'
        filename_code = 'capolar_k_116_n_128_ul.mat';
        S = open(filename_code);
        G = S.G;
        H = S.H;
        [k,n] = size(G);
        code_label = 'CAPOLAR';

    case 'CRC'
        k = 240;
        hex_poly = '0xd175';
        poly = koopman2matlab(hex_poly);
        [G,H,n] = make_CRC_GH(k,poly);
        code_label = ['CRC_' hex_poly];

    otherwise
        error('Unsupported code_class. Use ''CAPOLAR'' or ''CRC''.');
end

R = k/n;

%% Convert Eb/N0 to SNR
snr_db = ebn0 + 10*log10(R) + 10*log10(nmodbits);
numSNR = length(snr_db);

%% Pad modulation if ever needed
mod_pad = zeros(1,mod(n,nmodbits));

%% Determine number of bitonic compare-exchange steps
n_effective = 2^nextpow2(n);
m_bits      = log2(n_effective);
n_steps     = m_bits * (m_bits + 1) / 2;

%% Preallocate accumulators
tie_sum   = zeros(n_steps, numSNR);
total_sum = zeros(n_steps, numSNR);

% New: empirical input histogram of MSB3 values 0..7
msb_hist   = zeros(8, numSNR);
sample_sum = zeros(1, numSNR);

%% Monte Carlo loop
for ii = 1:numSNR

    sigma2 = 1/(10^(0.1*snr_db(ii)));

    awgnchan = comm.AWGNChannel( ...
        'NoiseMethod','Variance', ...
        'Variance',sigma2);

    fprintf('Running %s (n=%d, k=%d), Eb/N0 = %.2f dB, frames = %d\n', ...
        code_class, n, k, ebn0(ii), nFrames);

    for frm = 1:nFrames

        % Encode random information word
        u = binornd(ones(1,k),0.5);
        c = mod(u*G,2);

        % BPSK modulation + AWGN
        modOut = nrSymbolModulate([c mod_pad]', modulation);
        rSig   = awgnchan(modOut);

        % Soft demodulation
        y_soft = nrSymbolDemodulate(rSig, modulation, sigma2);
        y_soft = y_soft(1:n)';

        % ---------------------------------------------
        % Match the MEX quantization path exactly
        % ---------------------------------------------
        L = max(min(y_soft, 31), -31);        % clip to [-31,31]
        mag_q = round(abs(L));                % 5-bit magnitude 0...31
        msb3 = bitshift(uint8(mag_q), -2);    % keep top 3 MSBs -> 0...7

        for b = 0:7
            msb_hist(b+1,ii) = msb_hist(b+1,ii) + sum(msb3 == b);
        end
        sample_sum(ii) = sample_sum(ii) + numel(msb3);

        % Call precision-aware sorter MEX
        [~, tie_count, total_count] = pa_step_mex(y_soft);

        % Accumulate step statistics
        tie_sum(:,ii)   = tie_sum(:,ii)   + double(tie_count(:));
        total_sum(:,ii) = total_sum(:,ii) + double(total_count(:));
    end
end

%% Compute step-wise tie probabilities
tie_prob = tie_sum ./ total_sum;

% New: empirical MSB3 probability distribution
p_msb3_mc = msb_hist ./ sample_sum;

%% Store results in struct
result = struct();
result.class       = code_class;
result.label       = code_label;
result.n           = n;
result.k           = k;
result.R           = R;
result.G           = G;
result.H           = H;
result.ebn0        = ebn0;
result.snr_db      = snr_db;
result.nFrames     = nFrames;
result.n_effective = n_effective;
result.n_steps     = n_steps;
result.tie_sum     = tie_sum;
result.total_sum   = total_sum;
result.tie_prob    = tie_prob;
result.msb_hist    = msb_hist;
result.sample_sum  = sample_sum;
result.p_msb3_mc   = p_msb3_mc;

%% Save results
if strcmpi(code_class,'CRC')
    outname = sprintf('../ORBGRAND_Code/TIE_STATS_STEPS_%s_%d_%d.mat', code_label, n, k);
else
    outname = sprintf('../ORBGRAND_Code/TIE_STATS_STEPS_%s_%d_%d.mat', code_class, n, k);
end

save(outname, 'result');
fprintf('Saved: %s\n', outname);

%% Show empirical MSB3 distribution
disp('Empirical MATLAB p_msb3_mc (rows = bins 0..7, cols = Eb/N0 points):');
disp(p_msb3_mc);

%% Export for Python
writematrix(result.ebn0, 'ebn0.csv');
writematrix(result.tie_prob(1,:), 'tie_step1.csv');   % exact first CAE step
writematrix(result.R, 'rate.csv');
writematrix(result.p_msb3_mc, 'p_msb3_mc.csv');
fprintf('Exported: ebn0.csv, tie_step1.csv, rate.csv, p_msb3_mc.csv\n');

%% Quick plot: tie probability vs raw step
figure; clf;
hold on;
for ii = 1:numSNR
    plot(1:n_steps, result.tie_prob(:,ii), '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('$E_b/N_0 = %d$ dB', ebn0(ii)));
end
hold off;
grid on;
box on;
xlabel('Compare-exchange step','Interpreter','latex');
ylabel('Tie probability','Interpreter','latex');
title(sprintf('%s (%d,%d)', code_class, n, k), 'Interpreter','none');
legend('show','Interpreter','latex','Location','best');

%% Quick plot: tie probability vs normalized step
figure; clf;
xnorm = (1:n_steps) / n_steps;
hold on;
for ii = 1:numSNR
    plot(xnorm, result.tie_prob(:,ii), '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('$E_b/N_0 = %d$ dB', ebn0(ii)));
end
hold off;
grid on;
box on;
xlabel('Normalized compare-exchange step','Interpreter','latex');
ylabel('Tie probability','Interpreter','latex');
title(sprintf('%s (%d,%d)', code_class, n, k), 'Interpreter','none');
legend('show','Interpreter','latex','Location','best');