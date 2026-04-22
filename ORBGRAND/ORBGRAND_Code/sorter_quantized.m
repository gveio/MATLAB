% ==========================================================
% Monte Carlo framework for displacement statistics of
% quantized / reduced-precision / precision-aware sorters
%
% MEX required:
%   sorted_idx = sorter_quantized_mex(y_soft, mode)
%   sorted_idx = sorter_quantized_mex(y_soft, mode, n_sort_effective_forced, n_tie_effective_forced)
%
% mode mapping:
%   0 -> true-valued reference
%   1 -> baseline quantized
%   2 -> MSB2
%   3 -> MSB3
%   4 -> MSB4
%   5 -> PA1
%   6 -> PA2
%   7 -> PA1_pruned
%   8 -> PA2_pruned
%
% Output:
%   result.disp_pct
%   result.topk_overlap_pct
%   result.pa_compare
%
% Notes:
%   - Reference ordering is the true-valued abs(LLR) sorting.
%   - Displacement thresholds are evaluated ONLY over the
%     decoder-relevant top-LWmax elements of the reference ordering.
%   - Top-K overlap measures set agreement within that region.
%   - If force_hw_mapping = true, modes are evaluated using
%     sorter_quantized_mex(..., mode, 256, 256)
%     to match the hardware-oriented fixed 256-stage implementation.
% ==========================================================

clear; clc;

%% ----------------------------------------------------------
% User options
% -----------------------------------------------------------
code_class        = 'CRC';   % 'CAPOLAR' or 'CRC'
ebn0              = 7;       % can be scalar or vector
nFrames           = 2e4;     % Monte Carlo frames per SNR
LWmax             = 104;     % decoder-relevant top-K region
thresholds        = [0 1 2 3 5 10 20];

% Validate hardware-oriented implementation?
force_hw_mapping  = true;    % true -> use sorter_quantized_mex(y_soft, mode, 256, 256)
n_sort_forced     = 256;
n_tie_forced      = 256;

%% ----------------------------------------------------------
% Modulation: BPSK only
% -----------------------------------------------------------
modulation = 'BPSK';
nmodbits   = 1;

%% ----------------------------------------------------------
% Load / build selected code
% -----------------------------------------------------------
switch upper(code_class)

    case 'CAPOLAR'
        % CA-Polar(128,116)
        filename_code = 'capolar_k_116_n_128_ul.mat';
        S = open(filename_code);
        G = S.G;
        H = S.H;
        [k,n] = size(G);
        code_label = 'CAPOLAR';

    case 'CRC'
        % CRC(256,240)
        k = 240;
        hex_poly = '0xd175';
        poly = koopman2matlab(hex_poly);
        [G,H,n] = make_CRC_GH(k,poly);
        code_label = ['CRC_' hex_poly];

    otherwise
        error('Unsupported code_class. Use ''CAPOLAR'' or ''CRC''.');
end

R = k/n;

%% ----------------------------------------------------------
% Convert Eb/N0 to SNR
% -----------------------------------------------------------
snr_db = ebn0 + 10*log10(R) + 10*log10(nmodbits);
numSNR = length(snr_db);

%% ----------------------------------------------------------
% Padding for modulation if ever needed
% -----------------------------------------------------------
mod_pad = zeros(1,mod(n,nmodbits));

%% ----------------------------------------------------------
% Methods
% -----------------------------------------------------------
method_names = { ...
    'quant','msb2','msb3','msb4', ...
    'pa1','pa2', ...
    'pa1_pruned','pa2_pruned'};

method_modes = [1 2 3 4 5 6 7 8];
numMethods   = numel(method_names);

k_top  = min(LWmax, n);
numThr = numel(thresholds);

%% ----------------------------------------------------------
% Preallocate accumulators
% disp_count(threshold, method, snr) over TOP-K only
% topk_match(method, snr)
% -----------------------------------------------------------
disp_count  = zeros(numThr, numMethods, numSNR);
topk_match  = zeros(numMethods, numSNR);
total_elems = zeros(1, numSNR);   % counts k_top per frame

%% ----------------------------------------------------------
% Monte Carlo loop
% -----------------------------------------------------------
for ii = 1:numSNR

    sigma2 = 1/(10^(0.1*snr_db(ii)));

    awgnchan = comm.AWGNChannel( ...
        'NoiseMethod','Variance', ...
        'Variance',sigma2);

    fprintf('Running %s (n=%d, k=%d), Eb/N0 = %.2f dB, frames = %d\n', ...
        code_class, n, k, ebn0(ii), nFrames);

    for frm = 1:nFrames

        % --------------------------------------------------
        % Encode random information word
        % --------------------------------------------------
        u = binornd(ones(1,k),0.5);
        c = mod(u*G,2);

        % --------------------------------------------------
        % BPSK modulation + AWGN
        % --------------------------------------------------
        modOut = nrSymbolModulate([c mod_pad]', modulation);
        rSig   = awgnchan(modOut);

        % --------------------------------------------------
        % Soft demodulation
        % --------------------------------------------------
        y_soft = nrSymbolDemodulate(rSig, modulation, sigma2);
        y_soft = y_soft(1:n)';

        % --------------------------------------------------
        % Reference true-valued ordering
        % --------------------------------------------------
        idx_ref = double(sorter_quantized_mex(y_soft, 0));  % mode 0 = true-valued

        % Position map of reference ordering
        pos_ref = zeros(1,n);
        pos_ref(idx_ref) = 1:n;

        % Reference top-K set
        top_ref = idx_ref(1:k_top);

        % --------------------------------------------------
        % Compare all quantized / PA variants
        % --------------------------------------------------
        for mm = 1:numMethods

            if force_hw_mapping
                idx_test = double(sorter_quantized_mex( ...
                    y_soft, method_modes(mm), n_sort_forced, n_tie_forced));
            else
                idx_test = double(sorter_quantized_mex(y_soft, method_modes(mm)));
            end

            % Position map
            pos_test = zeros(1,n);
            pos_test(idx_test) = 1:n;

            % Element-wise displacement ONLY for the reference top-K elements
            disp_vec = abs(pos_test(top_ref) - pos_ref(top_ref));

            % Cumulative displacement counts over top-K only
            for tt = 1:numThr
                disp_count(tt,mm,ii) = disp_count(tt,mm,ii) + sum(disp_vec <= thresholds(tt));
            end

            % Top-K overlap
            top_test = idx_test(1:k_top);
            topk_match(mm,ii) = topk_match(mm,ii) + numel(intersect(top_ref, top_test));
        end

        total_elems(ii) = total_elems(ii) + k_top;
    end
end

%% ----------------------------------------------------------
% Normalize results
% -----------------------------------------------------------
disp_pct = zeros(numThr, numMethods, numSNR);
topk_overlap_pct = zeros(numMethods, numSNR);

for ii = 1:numSNR
    disp_pct(:,:,ii) = 100 * disp_count(:,:,ii) / total_elems(ii);
    topk_overlap_pct(:,ii) = 100 * topk_match(:,ii) / (nFrames * k_top);
end

%% ----------------------------------------------------------
% Direct PA vs PRUNED comparison
% -----------------------------------------------------------
result = struct();

pairs = {
    'pa1','pa1_pruned';
    'pa2','pa2_pruned'
};

pa_compare = struct();

for pp = 1:size(pairs,1)

    name_a = pairs{pp,1};
    name_b = pairs{pp,2};

    idx_a = find(strcmp(method_names, name_a));
    idx_b = find(strcmp(method_names, name_b));

    disp_diff = zeros(numThr, numSNR);
    topk_diff = zeros(1, numSNR);

    for ii = 1:numSNR
        disp_diff(:,ii) = disp_pct(:,idx_b,ii) - disp_pct(:,idx_a,ii);
        topk_diff(ii)   = topk_overlap_pct(idx_b,ii) - topk_overlap_pct(idx_a,ii);
    end

    pa_compare.(name_b).disp_diff = disp_diff;
    pa_compare.(name_b).topk_diff = topk_diff;
end

%% ----------------------------------------------------------
% Store results in struct
% -----------------------------------------------------------
result.class            = code_class;
result.label            = code_label;
result.n                = n;
result.k                = k;
result.R                = R;
result.G                = G;
result.H                = H;
result.ebn0             = ebn0;
result.snr_db           = snr_db;
result.nFrames          = nFrames;
result.thresholds       = thresholds;
result.k_top            = k_top;
result.method_names     = method_names;
result.method_modes     = method_modes;
result.disp_count       = disp_count;
result.disp_pct         = disp_pct;
result.topk_match       = topk_match;
result.topk_overlap_pct = topk_overlap_pct;
result.pa_compare       = pa_compare;
result.force_hw_mapping = force_hw_mapping;
result.n_sort_forced    = n_sort_forced;
result.n_tie_forced     = n_tie_forced;

%% ----------------------------------------------------------
% Save results
% -----------------------------------------------------------
if strcmpi(code_class,'CRC')
    outname = sprintf('../ORBGRAND_Code/DISP_STATS_%s_%d_%d.mat', code_label, n, k);
else
    outname = sprintf('../ORBGRAND_Code/DISP_STATS_%s_%d_%d.mat', code_class, n, k);
end

save(outname, 'result');
fprintf('Saved: %s\n', outname);

%% ----------------------------------------------------------
% Example printout for one representative SNR
% -----------------------------------------------------------
[~, idx_mid] = min(abs(ebn0 - 7));   % try to report around 7 dB

fprintf('\n====================================================\n');
fprintf('Displacement statistics at Eb/N0 = %.2f dB\n', ebn0(idx_mid));
fprintf('Reference = true-valued abs(LLR) ordering\n');
fprintf('Displacement evaluated over the reference Top-%d elements only\n', k_top);
fprintf('Top-K = %d\n', k_top);
if force_hw_mapping
    fprintf('Sorter evaluation = forced hardware mapping (%d,%d)\n', n_sort_forced, n_tie_forced);
else
    fprintf('Sorter evaluation = native mode\n');
end
fprintf('====================================================\n');

fprintf('%12s', 'Method');
for tt = 1:numThr
    if thresholds(tt) == 0
        fprintf('%10s', '=0');
    else
        fprintf('%10s', ['<=' num2str(thresholds(tt))]);
    end
end
fprintf('%12s\n', sprintf('Top-%d', k_top));

for mm = 1:numMethods
    fprintf('%12s', method_names{mm});
    for tt = 1:numThr
        fprintf('%10.2f', disp_pct(tt,mm,idx_mid));
    end
    fprintf('%12.2f\n', topk_overlap_pct(mm,idx_mid));
end

%% ----------------------------------------------------------
% Print pruned vs original deltas
% -----------------------------------------------------------
fprintf('\n====================================================\n');
fprintf('PRUNED vs ORIGINAL (delta)\n');
fprintf('====================================================\n');

for pp = 1:size(pairs,1)

    name_a = pairs{pp,1};
    name_b = pairs{pp,2};

    fprintf('\n%s vs %s at Eb/N0 = %.2f dB\n', name_b, name_a, ebn0(idx_mid));
    fprintf('Top-%d overlap delta: %.6f %%\n', ...
        k_top, result.pa_compare.(name_b).topk_diff(idx_mid));
    fprintf('=0 delta: %.6f %%\n', ...
        result.pa_compare.(name_b).disp_diff(1,idx_mid));
    fprintf('<=1 delta: %.6f %%\n', ...
        result.pa_compare.(name_b).disp_diff(2,idx_mid));
    fprintf('<=2 delta: %.6f %%\n', ...
        result.pa_compare.(name_b).disp_diff(3,idx_mid));
    fprintf('<=5 delta: %.6f %%\n', ...
        result.pa_compare.(name_b).disp_diff(5,idx_mid));
    fprintf('<=10 delta: %.6f %%\n', ...
        result.pa_compare.(name_b).disp_diff(6,idx_mid));
end

%% ----------------------------------------------------------
% Optional CSV export for the representative SNR
% -----------------------------------------------------------
table_mat = zeros(numMethods, numThr+1);
table_mat(:,1:numThr) = disp_pct(:,:,idx_mid).';
table_mat(:,numThr+1) = topk_overlap_pct(:,idx_mid);

csv_name = sprintf('disp_topk_table_%s_n%d_EbN0_%g.csv', code_class, n, ebn0(idx_mid));
writematrix(table_mat, csv_name);
fprintf('Saved CSV: %s\n', csv_name);

%% ----------------------------------------------------------
% Optional CSV export for PA vs PRUNED deltas
% -----------------------------------------------------------
delta_mat = zeros(size(pairs,1), numThr+1);
for pp = 1:size(pairs,1)
    name_b = pairs{pp,2};
    delta_mat(pp,1:numThr) = result.pa_compare.(name_b).disp_diff(:,idx_mid).';
    delta_mat(pp,numThr+1) = result.pa_compare.(name_b).topk_diff(idx_mid);
end

csv_name_delta = sprintf('disp_topk_delta_table_%s_n%d_EbN0_%g.csv', code_class, n, ebn0(idx_mid));
writematrix(delta_mat, csv_name_delta);
fprintf('Saved CSV: %s\n', csv_name_delta);

%% ----------------------------------------------------------
% Quick plot: top-K overlap vs Eb/N0
% -----------------------------------------------------------
figure; clf;
hold on;
for mm = 1:numMethods
    plot(ebn0, topk_overlap_pct(mm,:), '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', method_names{mm});
end
hold off;
grid on;
box on;
xlabel('$E_b/N_0$ (dB)', 'Interpreter','latex');
ylabel(sprintf('Top-%d overlap (%%)', k_top), 'Interpreter','latex');
title(sprintf('%s (%d,%d)', code_class, n, k), 'Interpreter','none');
legend('show', 'Interpreter','none', 'Location','best');

%% ----------------------------------------------------------
% Quick plot: exact-position match over Top-K only
% -----------------------------------------------------------
figure; clf;
hold on;
for mm = 1:numMethods
    plot(ebn0, squeeze(disp_pct(1,mm,:)), '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', method_names{mm});
end
hold off;
grid on;
box on;
xlabel('$E_b/N_0$ (dB)', 'Interpreter','latex');
ylabel(sprintf('Top-%d exact match (%%)', k_top), 'Interpreter','latex');
title(sprintf('%s (%d,%d)', code_class, n, k), 'Interpreter','none');
legend('show', 'Interpreter','none', 'Location','best');