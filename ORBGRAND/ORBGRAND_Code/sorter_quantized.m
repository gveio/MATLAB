% ==========================================================
% Monte Carlo framework for displacement statistics of
% quantized / reduced-precision / precision-aware sorters
%
% MEX required:
%   sorted_idx = sorter_quantized_mex(y_soft, mode)
%
% mode mapping:
%   0 -> true-valued reference
%   1 -> baseline quantized
%   2 -> MSB2
%   3 -> MSB3
%   4 -> MSB4
%   5 -> PA1
%   6 -> PA2
%
% Output:
%   result.disp_pct.(method)    : cumulative displacement %
%   result.topk_overlap.(method): top-LWmax overlap %
%   result.thresholds           : displacement thresholds
%
% Notes:
%   - Reference ordering is the true-valued abs(LLR) sorting.
%   - Quantized baseline shows quantization-induced deviation.
%   - Reduced-precision and PA variants show additional deviation.
% ==========================================================

clear; clc;

%% ----------------------------------------------------------
% User options
% -----------------------------------------------------------
code_class  = 'CRC';       % 'CAPOLAR' or 'CRC'
ebn0        = 7;       % Eb/N0 sweep
nFrames     = 2e4;         % Monte Carlo frames per SNR
LWmax       = 104;         % decoder-relevant top-K region
thresholds  = [0 1 2 3 5 10 20];

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
method_names = {'quant','msb2','msb3','msb4','pa1','pa2'};
method_modes = [1 2 3 4 5 6];
numMethods   = numel(method_names);

k_top = min(LWmax, n);
numThr = numel(thresholds);

%% ----------------------------------------------------------
% Preallocate accumulators
%
% disp_count(threshold, method, snr)
% topk_match(method, snr)
% -----------------------------------------------------------
disp_count = zeros(numThr, numMethods, numSNR);
topk_match = zeros(numMethods, numSNR);
total_elems = zeros(1, numSNR);

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

        % Position map of reference ordering:
        % pos_ref(bit_index) = sorted position
        pos_ref = zeros(1,n);
        pos_ref(idx_ref) = 1:n;

        % Top-K reference set
        top_ref = idx_ref(1:k_top);

        % --------------------------------------------------
        % Compare all quantized / PA variants
        % --------------------------------------------------
        for mm = 1:numMethods

            idx_test = double(sorter_quantized_mex(y_soft, method_modes(mm)));

            % Position map
            pos_test = zeros(1,n);
            pos_test(idx_test) = 1:n;

            % Element-wise displacement
            disp_vec = abs(pos_test - pos_ref);

            % Cumulative displacement counts
            for tt = 1:numThr
                disp_count(tt,mm,ii) = disp_count(tt,mm,ii) + sum(disp_vec <= thresholds(tt));
            end

            % Top-K overlap
            top_test = idx_test(1:k_top);
            topk_match(mm,ii) = topk_match(mm,ii) + numel(intersect(top_ref, top_test));
        end

        total_elems(ii) = total_elems(ii) + n;
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
% Store results in struct
% -----------------------------------------------------------
result = struct();
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
fprintf('Top-K = %d\n', k_top);
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
% Optional CSV export for the representative SNR
% -----------------------------------------------------------
table_mat = zeros(numMethods, numThr+1);
table_mat(:,1:numThr) = disp_pct(:,:,idx_mid).';
table_mat(:,numThr+1) = topk_overlap_pct(:,idx_mid);

csv_name = sprintf('disp_table_%s_n%d_EbN0_%g.csv', code_class, n, ebn0(idx_mid));
writematrix(table_mat, csv_name);
fprintf('Saved CSV: %s\n', csv_name);

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
% Quick plot: exact-position match (displacement = 0)
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
ylabel('Exact position match (\%)', 'Interpreter','latex');
title(sprintf('%s (%d,%d)', code_class, n, k), 'Interpreter','none');
legend('show', 'Interpreter','none', 'Location','best');