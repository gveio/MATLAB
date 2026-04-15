% ==========================================================
% Monte Carlo framework for native-128 vs mapped-to-256
% sorter comparison
%
% Requires updated MEX:
%   sorted_idx = sorter_quantized_mex(y_soft, mode)
%   sorted_idx = sorter_quantized_mex(y_soft, mode, 256)
%
% Modes:
%   1 -> baseline quantized
%   2 -> MSB2
%   3 -> MSB3
%   4 -> MSB4
%   5 -> PA1
%   6 -> PA2
%
% Reference:
%   native 128 sorter output for the same mode
%
% Test:
%   128 inputs mapped onto 256 sorting network
% ==========================================================

clear; clc;

%% ----------------------------------------------------------
% User options
% -----------------------------------------------------------
code_class  = 'CAPOLAR';   % this check is mainly for n=128
ebn0        = 7;           % representative point
nFrames     = 2e4;
LWmax       = 104;
thresholds  = [0 1 2 3 5 10 20];

%% ----------------------------------------------------------
% Modulation
% -----------------------------------------------------------
modulation = 'BPSK';
nmodbits   = 1;

%% ----------------------------------------------------------
% Load / build code
% -----------------------------------------------------------
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
        error('Unsupported code_class.');
end

if n ~= 128
    warning('This script is intended for n=128 vs mapped-256 validation.');
end

R = k/n;
snr_db = ebn0 + 10*log10(R) + 10*log10(nmodbits);
sigma2 = 1/(10^(0.1*snr_db));

%% ----------------------------------------------------------
% Padding for modulation if needed
% -----------------------------------------------------------
mod_pad = zeros(1,mod(n,nmodbits));

%% ----------------------------------------------------------
% Modes to compare
% -----------------------------------------------------------
method_names = {'baseline','msb2','msb3','msb4','pa1','pa2'};
method_modes = [1 2 3 4 5 6];
numMethods   = numel(method_names);

k_top = min(LWmax, n);
numThr = numel(thresholds);

%% ----------------------------------------------------------
% Accumulators
% -----------------------------------------------------------
disp_count = zeros(numThr, numMethods);
topk_match = zeros(numMethods, 1);
total_elems = 0;

%% ----------------------------------------------------------
% AWGN channel
% -----------------------------------------------------------
awgnchan = comm.AWGNChannel( ...
    'NoiseMethod','Variance', ...
    'Variance',sigma2);

fprintf('Running native-vs-mapped validation for %s (n=%d, k=%d), Eb/N0 = %.2f dB, frames = %d\n', ...
    code_class, n, k, ebn0, nFrames);

%% ----------------------------------------------------------
% Monte Carlo loop
% -----------------------------------------------------------
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

    for mm = 1:numMethods
        mode = method_modes(mm);

        % Native 128 sorter output (reference)
        idx_native = double(sorter_quantized_mex(y_soft, mode));

        % 128 mapped onto 256 sorting network
        idx_map256 = double(sorter_quantized_mex(y_soft, mode, 256));

        % Position maps
        pos_native = zeros(1,n);
        pos_native(idx_native) = 1:n;

        pos_map256 = zeros(1,n);
        pos_map256(idx_map256) = 1:n;

        % Displacement relative to native implementation
        disp_vec = abs(pos_map256 - pos_native);

        % Threshold counts
        for tt = 1:numThr
            disp_count(tt,mm) = disp_count(tt,mm) + sum(disp_vec <= thresholds(tt));
        end

        % Top-LWmax overlap
        top_native = idx_native(1:k_top);
        top_map256 = idx_map256(1:k_top);
        topk_match(mm) = topk_match(mm) + numel(intersect(top_native, top_map256));
    end

    total_elems = total_elems + n;
end

%% ----------------------------------------------------------
% Normalize
% -----------------------------------------------------------
disp_pct = 100 * disp_count / total_elems;
topk_overlap_pct = 100 * topk_match / (nFrames * k_top);

%% ----------------------------------------------------------
% Print results
% -----------------------------------------------------------
fprintf('\n====================================================\n');
fprintf('Native 128 vs mapped-to-256 displacement statistics\n');
fprintf('Reference = native 128 sorter output (same mode)\n');
fprintf('Top-K = %d\n', k_top);
fprintf('====================================================\n');

fprintf('%12s', 'Method');
for i = 1:numel(thresholds)
    if thresholds(i) == 0
        fprintf('%10s', '=0');
    else
        fprintf('%10s', ['<=' num2str(thresholds(i))]);
    end
end
fprintf('%12s\n', sprintf('Top-%d', k_top));

for mm = 1:numMethods
    fprintf('%12s', method_names{mm});
    fprintf('%10.2f', disp_pct(:,mm));
    fprintf('%12.2f\n', topk_overlap_pct(mm));
end

%% ----------------------------------------------------------
% Store results
% -----------------------------------------------------------
result = struct();
result.class            = code_class;
result.label            = code_label;
result.n                = n;
result.k                = k;
result.R                = R;
result.ebn0             = ebn0;
result.snr_db           = snr_db;
result.nFrames          = nFrames;
result.thresholds       = thresholds;
result.k_top            = k_top;
result.method_names     = method_names;
result.method_modes     = method_modes;
result.disp_count       = disp_count;
result.disp_pct         = disp_pct;
result.topk_overlap_pct = topk_overlap_pct;

outname = sprintf('MAP256_VALIDATION_%s_%d_%d_EbN0_%g.mat', code_label, n, k, ebn0);
save(outname, 'result');
fprintf('Saved: %s\n', outname);

%% ----------------------------------------------------------
% Optional CSV
% -----------------------------------------------------------
table_mat = [disp_pct.' topk_overlap_pct];
csv_name = sprintf('map256_validation_%s_n%d_EbN0_%g.csv', code_class, n, ebn0);
writematrix(table_mat, csv_name);
fprintf('Saved CSV: %s\n', csv_name);