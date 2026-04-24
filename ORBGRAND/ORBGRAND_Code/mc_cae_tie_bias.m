% ==========================================================
% Monte Carlo per-CAE MSB tie + tie-swap bias analysis
%
% MEX required:
%   [sorted_idx, tie_count, total_count, tie_swap_count] = ...
%       pa_cae_tie_bias_mex(y_soft)
%
% Outputs:
%   result.p_swap_given_tie(step, cae, snr)
%   result.bias_from_half(step, cae, snr)
% ==========================================================

clear; clc;

%% User options
code_class = 'CRC';     % 'CAPOLAR' or 'CRC'
ebn0       = 4:1:7;
nFrames    = 2e4;

modulation = 'BPSK';
nmodbits   = 1;

%% Load code
switch upper(code_class)

    case 'CAPOLAR'
        S = open('capolar_k_116_n_128_ul.mat');
        G = S.G; H = S.H;
        [k,n] = size(G);
        code_label = 'CAPOLAR';

    case 'CRC'
        k = 240;
        hex_poly = '0xd175';
        poly = koopman2matlab(hex_poly);
        [G,H,n] = make_CRC_GH(k,poly);
        code_label = ['CRC_' hex_poly];

    otherwise
        error('Unsupported code_class');
end

R = k/n;

%% SNR
snr_db = ebn0 + 10*log10(R) + 10*log10(nmodbits);
numSNR = length(snr_db);

mod_pad = zeros(1,mod(n,nmodbits));

%% Bitonic dimensions
n_effective = 2^nextpow2(n);
m_bits      = log2(n_effective);
n_steps     = m_bits*(m_bits+1)/2;
n_cae       = n_effective/2;

%% Accumulators
tie_sum        = zeros(n_steps, n_cae, numSNR);
total_sum      = zeros(n_steps, n_cae, numSNR);
tie_swap_sum   = zeros(n_steps, n_cae, numSNR);

%% Monte Carlo
for ii = 1:numSNR

    sigma2 = 1/(10^(0.1*snr_db(ii)));
    awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);

    fprintf('Running %s Eb/N0=%.2f dB\n', code_class, ebn0(ii));

    for frm = 1:nFrames

        u = binornd(ones(1,k),0.5);
        c = mod(u*G,2);

        modOut = nrSymbolModulate([c mod_pad]', modulation);
        rSig   = awgnchan(modOut);

        y_soft = nrSymbolDemodulate(rSig, modulation, sigma2);
        y_soft = y_soft(1:n)';

        [~, tie_c, total_c, tie_swap_c] = pa_cae_tie_bias_mex(y_soft);

        tie_sum(:,:,ii)      = tie_sum(:,:,ii)      + double(tie_c);
        total_sum(:,:,ii)    = total_sum(:,:,ii)    + double(total_c);
        tie_swap_sum(:,:,ii) = tie_swap_sum(:,:,ii) + double(tie_swap_c);
    end
end

%% Compute probabilities
p_swap_given_tie = nan(size(tie_sum));
valid = tie_sum > 0;
p_swap_given_tie(valid) = tie_swap_sum(valid) ./ tie_sum(valid);

bias_from_half = abs(p_swap_given_tie - 0.5);

%% Store
result = struct();
result.p_swap_given_tie = p_swap_given_tie;
result.bias_from_half   = bias_from_half;
result.tie_sum          = tie_sum;
result.total_sum        = total_sum;
result.tie_swap_sum     = tie_swap_sum;
result.ebn0             = ebn0;
result.n_steps          = n_steps;
result.n_cae            = n_cae;

%% Save
save('CAE_TIE_BIAS.mat','result');

%% ==========================================================
%% PLOTTING
%% ==========================================================

ii = 4; % choose SNR index (e.g. 6 dB)

%% 1) Swap probability heatmap
figure; clf;
imagesc(p_swap_given_tie(:,:,ii));
colorbar;
caxis([0 1]);
xlabel('CAE index','Interpreter','latex');
ylabel('Bitonic step','Interpreter','latex');
title(sprintf('$P(\\mathrm{swap}|\\mathrm{MSB\\ tie}), E_b/N_0=%d$ dB', ebn0(ii)), ...
    'Interpreter','latex');

%% 2) Bias magnitude heatmap
figure; clf;
imagesc(bias_from_half(:,:,ii));
colorbar;
caxis([0 0.5]);
xlabel('CAE index','Interpreter','latex');
ylabel('Bitonic step','Interpreter','latex');
title(sprintf('Bias magnitude, $E_b/N_0=%d$ dB', ebn0(ii)), ...
    'Interpreter','latex');

%% 3) Histogram of bias values (optional insight)
figure; clf;
histogram(p_swap_given_tie(:,:,ii), 20);
xlabel('$P(\mathrm{swap}|\mathrm{tie})$','Interpreter','latex');
ylabel('Count');
title('Distribution across all CAEs','Interpreter','latex');
grid on;