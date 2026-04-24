% ==========================================================
% Monte Carlo framework for step-wise MSB tie probability
% and tie-swap bias using precision-aware bitonic sorter MEX
%
% MEX required:
%   [sorted_idx, tie_count, total_count, tie_swap_count] = ...
%       pa_step_tie_bias_mex(y_soft)
%
% Main metrics:
%   tie_prob(t) =
%       tie_sum(t) / total_sum(t)
%
%   p_swap_given_tie(t) =
%       tie_swap_sum(t) / tie_sum(t)
%
%   bias_from_half(t) =
%       abs(p_swap_given_tie(t) - 0.5)
% ==========================================================

clear; clc;

%% User options
code_class = 'CRC';     % 'CAPOLAR' or 'CRC'
ebn0       = 4:1:7;
nFrames    = 2e4;

%% Modulation
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

%% Padding
mod_pad = zeros(1,mod(n,nmodbits));

%% Number of bitonic steps
n_effective = 2^nextpow2(n);
m_bits      = log2(n_effective);
n_steps     = m_bits * (m_bits + 1) / 2;

%% Accumulators
tie_sum        = zeros(n_steps, numSNR);
total_sum      = zeros(n_steps, numSNR);
tie_swap_sum   = zeros(n_steps, numSNR);

%% Monte Carlo loop
for ii = 1:numSNR

    sigma2 = 1/(10^(0.1*snr_db(ii)));

    awgnchan = comm.AWGNChannel( ...
        'NoiseMethod','Variance', ...
        'Variance',sigma2);

    fprintf('Running %s (n=%d, k=%d), Eb/N0 = %.2f dB, frames = %d\n', ...
        code_class, n, k, ebn0(ii), nFrames);

    for frm = 1:nFrames

        u = binornd(ones(1,k),0.5);
        c = mod(u*G,2);

        modOut = nrSymbolModulate([c mod_pad]', modulation);
        rSig   = awgnchan(modOut);

        y_soft = nrSymbolDemodulate(rSig, modulation, sigma2);
        y_soft = y_soft(1:n)';

        [~, tie_count, total_count, tie_swap_count] = ...
            pa_step_tie_bias_mex(y_soft);

        tie_sum(:,ii)      = tie_sum(:,ii)      + double(tie_count(:));
        total_sum(:,ii)    = total_sum(:,ii)    + double(total_count(:));
        tie_swap_sum(:,ii) = tie_swap_sum(:,ii) + double(tie_swap_count(:));
    end
end

%% Metrics
tie_prob = tie_sum ./ total_sum;

p_swap_given_tie = nan(size(tie_sum));
valid = tie_sum > 0;
p_swap_given_tie(valid) = tie_swap_sum(valid) ./ tie_sum(valid);

bias_from_half = abs(p_swap_given_tie - 0.5);

%% Store results
result = struct();
result.class             = code_class;
result.label             = code_label;
result.n                 = n;
result.k                 = k;
result.R                 = R;
result.G                 = G;
result.H                 = H;
result.ebn0              = ebn0;
result.snr_db            = snr_db;
result.nFrames           = nFrames;
result.n_effective       = n_effective;
result.n_steps           = n_steps;
result.tie_sum           = tie_sum;
result.total_sum         = total_sum;
result.tie_swap_sum      = tie_swap_sum;
result.tie_prob          = tie_prob;
result.p_swap_given_tie  = p_swap_given_tie;
result.bias_from_half    = bias_from_half;

%% Save
if strcmpi(code_class,'CRC')
    outname = sprintf('../ORBGRAND_Code/TIE_BIAS_STEPS_%s_%d_%d.mat', ...
        code_label, n, k);
else
    outname = sprintf('../ORBGRAND_Code/TIE_BIAS_STEPS_%s_%d_%d.mat', ...
        code_class, n, k);
end

save(outname, 'result');
fprintf('Saved: %s\n', outname);

%% Export
writematrix(result.ebn0, 'tie_bias_ebn0.csv');
writematrix(result.tie_prob, 'tie_prob_steps.csv');
writematrix(result.p_swap_given_tie, 'p_swap_given_tie_steps.csv');
writematrix(result.bias_from_half, 'bias_from_half_steps.csv');
writematrix(result.R, 'tie_bias_rate.csv');

%% Plot 1: tie probability
figure; clf;
hold on;
for ii = 1:numSNR
    plot(1:n_steps, result.tie_prob(:,ii), '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('$E_b/N_0 = %d$ dB', ebn0(ii)));
end
hold off;
grid on; box on;
xlabel('Bitonic compare-exchange step','Interpreter','latex');
ylabel('$P(\mathrm{MSB\ tie})$','Interpreter','latex');
title(sprintf('%s (%d,%d): MSB tie probability', code_class, n, k), ...
    'Interpreter','none');
legend('show','Interpreter','latex','Location','best');

%% Plot 2: swap probability conditioned on MSB tie
figure; clf;
hold on;
for ii = 1:numSNR
    plot(1:n_steps, result.p_swap_given_tie(:,ii), '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('$E_b/N_0 = %d$ dB', ebn0(ii)));
end
yline(0.5, '--', 'LineWidth', 1.2, ...
    'DisplayName', '$0.5$');
hold off;
grid on; box on;
ylim([0 1]);
xlabel('Bitonic compare-exchange step','Interpreter','latex');
ylabel('$P(\mathrm{swap}\mid\mathrm{MSB\ tie})$','Interpreter','latex');
title(sprintf('%s (%d,%d): tie-swap bias', code_class, n, k), ...
    'Interpreter','none');
legend('show','Interpreter','latex','Location','best');

%% Plot 3: absolute bias from 0.5
figure; clf;
hold on;
for ii = 1:numSNR
    plot(1:n_steps, result.bias_from_half(:,ii), '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('$E_b/N_0 = %d$ dB', ebn0(ii)));
end
hold off;
grid on; box on;
ylim([0 0.5]);
xlabel('Bitonic compare-exchange step','Interpreter','latex');
ylabel('$|P(\mathrm{swap}\mid\mathrm{MSB\ tie}) - 0.5|$', ...
    'Interpreter','latex');
title(sprintf('%s (%d,%d): tie-decision bias magnitude', code_class, n, k), ...
    'Interpreter','none');
legend('show','Interpreter','latex','Location','best');