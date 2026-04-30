% ==========================================================
% AR(1) correlated-noise sweep for MSB3-tie full-precision behavior
%
% Measures:
%   P(full swap | MSB3 tie)
%   P(full keep | MSB3 tie)
%
% Also includes AWGN-style diagnostics for a selected rho:
%   1) Heatmap of P(full swap | MSB3 tie)
%   2) Step/distance probability summary
%   3) Last-stage distance contribution
%   4) CAE concentration curve
% ==========================================================

clear; clc;

%% User options
code_class = 'CRC';     % 'CAPOLAR' or 'CRC'
ebn0       = 7;
nFrames    = 3e4;

rho_list = [0 0.3 0.5 0.7 0.9 -0.3 -0.5 -0.7 -0.9];

rho_plot = 0.9;         % rho used for AWGN-style diagnostic plots

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
sigma2 = 1/(10^(0.1*snr_db));

mod_pad = zeros(1,mod(n,nmodbits));

%% Bitonic metadata
n_effective = 2^nextpow2(n);
m_bits      = log2(n_effective);
n_steps     = m_bits*(m_bits+1)/2;
n_cae       = n_effective/2;

step_outer = zeros(n_steps,1);
step_local = zeros(n_steps,1);
step_dist  = zeros(n_steps,1);

cnt = 0;
for s = 1:m_bits
    for kstep = 1:s
        cnt = cnt + 1;
        step_outer(cnt) = s;
        step_local(cnt) = kstep;
        step_dist(cnt)  = 2^(s-kstep);
    end
end

last_stage_steps = find(step_outer == m_bits);

last2_outer_start = (m_bits-2)*(m_bits-1)/2 + 1;
last2_steps       = last2_outer_start:n_steps;

last2_dist_labels = strings(numel(last2_steps),1);
for i = 1:numel(last2_steps)
    st = last2_steps(i);
    last2_dist_labels(i) = sprintf('step %d, d=%d', st, step_dist(st));
end

fprintf('n = %d, n_effective = %d, steps = %d\n', n, n_effective, n_steps);
fprintf('Eb/N0 = %.2f dB, SNR = %.4f dB, sigma2 = %.6g\n', ebn0, snr_db, sigma2);
fprintf('Last two stages: steps %d to %d\n', last2_steps(1), last2_steps(end));

%% Results
nrho = length(rho_list);

P_tie_region            = zeros(nrho,1);
P_swap_given_tie_region = zeros(nrho,1);
P_keep_given_tie_region = zeros(nrho,1);
P_tie_and_swap_region   = zeros(nrho,1);

P_swap_step = nan(n_steps,nrho);
P_tie_step  = zeros(n_steps,nrho);

%% Full per-CAE storage for AWGN-style diagnostics
tie_sum_all      = zeros(n_steps,n_cae,nrho);
total_sum_all    = zeros(n_steps,n_cae,nrho);
tie_swap_sum_all = zeros(n_steps,n_cae,nrho);

%% Main sweep
for rid = 1:nrho

    rho = rho_list(rid);

    fprintf('\n========================================\n');
    fprintf('Running rho = %.2f\n', rho);
    fprintf('========================================\n');

    tie_sum      = zeros(n_steps,n_cae);
    total_sum    = zeros(n_steps,n_cae);
    tie_swap_sum = zeros(n_steps,n_cae);

    for frm = 1:nFrames

        if mod(frm,1000) == 0
            fprintf('  Frame %d / %d\n', frm, nFrames);
        end

        %% Encode
        u = binornd(ones(1,k),0.5);
        c = mod(u*G,2);

        %% MATLAB BPSK modulation
        modOut = nrSymbolModulate([c mod_pad]', modulation);
        modOut = modOut(1:n).';

        %% AR(1)-correlated Gaussian noise
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

        %% MATLAB BPSK demodulation
        y_soft = nrSymbolDemodulate(rSig.', modulation, sigma2);
        y_soft = y_soft(1:n).';

        %% MEX diagnostic
        [~, tie_c, total_c, tie_swap_c] = pa_cae_tie_bias_mex(y_soft);

        tie_sum      = tie_sum      + double(tie_c);
        total_sum    = total_sum    + double(total_c);
        tie_swap_sum = tie_swap_sum + double(tie_swap_c);
    end

    %% Store full per-CAE counters
    tie_sum_all(:,:,rid)      = tie_sum;
    total_sum_all(:,:,rid)    = total_sum;
    tie_swap_sum_all(:,:,rid) = tie_swap_sum;

    %% Region-level summary: last two stages
    total_ties_region  = sum(tie_sum(last2_steps,:),'all');
    total_swaps_region = sum(tie_swap_sum(last2_steps,:),'all');
    total_comp_region  = sum(total_sum(last2_steps,:),'all');

    P_tie_region(rid)            = total_ties_region / total_comp_region;
    P_swap_given_tie_region(rid) = total_swaps_region / total_ties_region;
    P_keep_given_tie_region(rid) = 1 - P_swap_given_tie_region(rid);
    P_tie_and_swap_region(rid)   = total_swaps_region / total_comp_region;

    %% Step-level summary
    for st = 1:n_steps
        total_ties_st  = sum(tie_sum(st,:),'all');
        total_swaps_st = sum(tie_swap_sum(st,:),'all');
        total_comp_st  = sum(total_sum(st,:),'all');

        P_tie_step(st,rid) = total_ties_st / total_comp_st;

        if total_ties_st > 0
            P_swap_step(st,rid) = total_swaps_st / total_ties_st;
        end
    end

    fprintf('\nrho = %.2f summary, last two stages:\n', rho);
    fprintf('P(tie)                  = %.4f\n', P_tie_region(rid));
    fprintf('P(full swap | MSB tie)  = %.4f\n', P_swap_given_tie_region(rid));
    fprintf('P(full keep | MSB tie)  = %.4f\n', P_keep_given_tie_region(rid));
    fprintf('P(tie and full swap)    = %.4f\n', P_tie_and_swap_region(rid));
end

%% Summary table
Summary = table( ...
    rho_list(:), ...
    P_tie_region, ...
    P_swap_given_tie_region, ...
    P_keep_given_tie_region, ...
    P_tie_and_swap_region, ...
    'VariableNames', ...
    {'rho','P_tie','P_full_swap_given_MSB_tie', ...
     'P_full_keep_given_MSB_tie','P_tie_and_full_swap'} ...
);

disp(Summary);

save('AR1_MSB3_TIE_SWAP_SWEEP_WITH_DIAGNOSTICS.mat', ...
    'Summary', 'P_swap_step', 'P_tie_step', ...
    'tie_sum_all', 'total_sum_all', 'tie_swap_sum_all', ...
    'rho_list', 'step_outer', 'step_local', 'step_dist', ...
    'last2_steps', 'last_stage_steps', 'code_label', 'ebn0', ...
    'snr_db', 'sigma2');

%% ==========================================================
%% Sweep plots
%% ==========================================================

%% Plot 1: swap/keep ratio vs rho
figure; clf;
plot(rho_list, P_swap_given_tie_region, '-o', 'LineWidth', 1.5); hold on;
plot(rho_list, P_keep_given_tie_region, '-s', 'LineWidth', 1.5);
grid on;
xlabel('$\rho$','Interpreter','latex');
ylabel('Conditional probability','Interpreter','latex');
legend({'$P(\mathrm{full\ swap}\mid\mathrm{MSB3\ tie})$', ...
        '$P(\mathrm{full\ keep}\mid\mathrm{MSB3\ tie})$'}, ...
        'Interpreter','latex', 'Location','best');
title(sprintf('MSB3-tie resolution bias, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% Plot 2: tie probability vs rho
figure; clf;
plot(rho_list, P_tie_region, '-o', 'LineWidth', 1.5); hold on;
plot(rho_list, P_tie_and_swap_region, '-s', 'LineWidth', 1.5);
grid on;
xlabel('$\rho$','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
legend({'$P(\mathrm{MSB3\ tie})$', ...
        '$P(\mathrm{MSB3\ tie\ and\ full\ swap})$'}, ...
        'Interpreter','latex', 'Location','best');
title(sprintf('Tie occurrence under AR(1) correlated noise, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% Plot 3: step-level swap-given-tie in last two stages
figure; clf;
hold on;

for rid = 1:nrho
    plot(last2_steps, P_swap_step(last2_steps,rid), '-o', ...
        'DisplayName', sprintf('\\rho=%.1f', rho_list(rid)));
end

grid on;
xlabel('Bitonic step','Interpreter','latex');
ylabel('$P(\mathrm{full\ swap}\mid\mathrm{MSB3\ tie})$', ...
    'Interpreter','latex');
legend('Location','best');
title(sprintf('Step-level full-swap probability under MSB3 ties, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% ==========================================================
%% AWGN-style diagnostics for selected rho
%% ==========================================================

[~, rid_plot] = min(abs(rho_list - rho_plot));
rho_sel = rho_list(rid_plot);

tie_sum_r      = tie_sum_all(:,:,rid_plot);
total_sum_r    = total_sum_all(:,:,rid_plot);
tie_swap_sum_r = tie_swap_sum_all(:,:,rid_plot);

p_tie_r = zeros(size(tie_sum_r));
valid_total = total_sum_r > 0;
p_tie_r(valid_total) = tie_sum_r(valid_total) ./ total_sum_r(valid_total);

p_swap_given_tie_r = nan(size(tie_sum_r));
valid_tie = tie_sum_r > 0;
p_swap_given_tie_r(valid_tie) = tie_swap_sum_r(valid_tie) ./ tie_sum_r(valid_tie);

p_tie_and_swap_r = zeros(size(tie_sum_r));
p_tie_and_swap_r(valid_total) = tie_swap_sum_r(valid_total) ./ total_sum_r(valid_total);

%% Last-two-stage step summary
rows = [];
total_last2_swaps = sum(tie_swap_sum_r(last2_steps,:),'all');

for st = last2_steps

    total_ties_st  = sum(tie_sum_r(st,:),'all');
    total_swaps_st = sum(tie_swap_sum_r(st,:),'all');
    total_comp_st  = sum(total_sum_r(st,:),'all');

    P_tie_st = total_ties_st / total_comp_st;
    P_swap_given_tie_st = total_swaps_st / total_ties_st;
    P_tie_and_swap_st = total_swaps_st / total_comp_st;
    frac_last2_swaps_st = total_swaps_st / total_last2_swaps;

    rows = [rows; ...
        st, step_outer(st), step_local(st), step_dist(st), ...
        total_ties_st, total_swaps_st, total_comp_st, ...
        P_tie_st, P_swap_given_tie_st, P_tie_and_swap_st, ...
        frac_last2_swaps_st];
end

Last2Summary_r = array2table(rows, ...
    'VariableNames', ...
    {'step','outer_stage','local_step','distance', ...
     'total_ties','full_swaps','total_comparisons', ...
     'P_tie','P_full_swap_given_tie','P_tie_and_full_swap', ...
     'fraction_of_last2_full_swaps'});

fprintf('\nLAST-TWO-STAGE SUMMARY FOR rho = %.2f:\n', rho_sel);
disp(Last2Summary_r);

%% Last-stage distance summary
rows = [];
total_last_stage_swaps = sum(tie_swap_sum_r(last_stage_steps,:),'all');

for st = last_stage_steps(:)'

    d = step_dist(st);

    total_ties_d  = sum(tie_sum_r(st,:),'all');
    total_swaps_d = sum(tie_swap_sum_r(st,:),'all');
    total_comp_d  = sum(total_sum_r(st,:),'all');

    P_tie_d = total_ties_d / total_comp_d;
    P_swap_given_tie_d = total_swaps_d / total_ties_d;
    P_tie_and_swap_d = total_swaps_d / total_comp_d;
    frac_last_stage_swaps_d = total_swaps_d / total_last_stage_swaps;

    rows = [rows; ...
        d, st, total_ties_d, total_swaps_d, total_comp_d, ...
        P_tie_d, P_swap_given_tie_d, P_tie_and_swap_d, ...
        frac_last_stage_swaps_d];
end

LastStageSummary_r = array2table(rows, ...
    'VariableNames', ...
    {'distance','step','total_ties','full_swaps','total_comparisons', ...
     'P_tie','P_full_swap_given_tie','P_tie_and_full_swap', ...
     'fraction_of_last_stage_full_swaps'});

LastStageSummary_r = sortrows(LastStageSummary_r, 'distance', 'ascend');

fprintf('\nLAST-STAGE DISTANCE SUMMARY FOR rho = %.2f:\n', rho_sel);
disp(LastStageSummary_r);

%% Uniformity check
valid_last2_psgt = p_swap_given_tie_r(last2_steps,:);
valid_last2_psgt = valid_last2_psgt(~isnan(valid_last2_psgt));

mean_psgt = mean(valid_last2_psgt);
std_psgt  = std(valid_last2_psgt);
min_psgt  = min(valid_last2_psgt);
max_psgt  = max(valid_last2_psgt);

fprintf('\nUNIFORMITY OF P(full swap | MSB3 tie), rho = %.2f:\n', rho_sel);
fprintf('mean = %.4f\n', mean_psgt);
fprintf('std  = %.4f\n', std_psgt);
fprintf('min  = %.4f\n', min_psgt);
fprintf('max  = %.4f\n', max_psgt);

%% Correlation between tie probability and full-swap fraction
corr_last2 = corr(Last2Summary_r.P_tie, Last2Summary_r.fraction_of_last2_full_swaps);
corr_last_stage = corr(LastStageSummary_r.P_tie, LastStageSummary_r.fraction_of_last_stage_full_swaps);

fprintf('\nDISTANCE EFFECT CHECK, rho = %.2f:\n', rho_sel);
fprintf('corr(P_tie, full-swap fraction) last two stages = %.4f\n', corr_last2);
fprintf('corr(P_tie, full-swap fraction) last stage      = %.4f\n', corr_last_stage);

%% CAE concentration ranking
last2_mask = false(n_steps,n_cae);
last2_mask(last2_steps,:) = true;

valid = last2_mask & (tie_swap_sum_r > 0);
[step_id, cae_id] = find(valid);

T_last2_r = table( ...
    step_id, ...
    step_outer(step_id), ...
    step_local(step_id), ...
    step_dist(step_id), ...
    cae_id, ...
    tie_sum_r(valid), ...
    tie_swap_sum_r(valid), ...
    total_sum_r(valid), ...
    p_tie_r(valid), ...
    p_swap_given_tie_r(valid), ...
    p_tie_and_swap_r(valid), ...
    'VariableNames', ...
    {'step','outer_stage','local_step','distance','cae', ...
     'tie_count','full_swap_count','total_count', ...
     'P_tie','P_full_swap_given_tie','P_tie_and_full_swap'} ...
);

T_last2_r = sortrows(T_last2_r, 'full_swap_count', 'descend');

total_swaps_last2 = sum(tie_swap_sum_r(last2_steps,:),'all');
T_last2_r.cum_full_swaps = cumsum(T_last2_r.full_swap_count);
T_last2_r.cum_fraction = T_last2_r.cum_full_swaps / total_swaps_last2;

fprintf('\nCAE CONCENTRATION FOR rho = %.2f:\n', rho_sel);
fprintf('Total CAEs in last two stages : %d\n', numel(last2_steps)*n_cae);
fprintf('CAEs with full swaps          : %d\n', height(T_last2_r));
fprintf('Total full swaps              : %d\n', total_swaps_last2);

Klist = [10 50 100 200 400 800 1200 1600];

fprintf('\nTop-K cumulative full-swap fraction:\n');
for kk = 1:length(Klist)
    Ksel = Klist(kk);
    if Ksel <= height(T_last2_r)
        fprintf('Top %4d -> fraction=%.4f\n', Ksel, T_last2_r.cum_fraction(Ksel));
    end
end

%% Diagnostic Plot 4: heatmap
figure; clf;
imagesc(p_swap_given_tie_r(last2_steps,:));
colorbar;
caxis([0 1]);
xlabel('CAE index','Interpreter','latex');
ylabel('Last-two-stage step','Interpreter','latex');
yticks(1:numel(last2_steps));
yticklabels(last2_dist_labels);
title(sprintf('$P(\\mathrm{full\\ swap}\\mid\\mathrm{MSB3\\ tie})$, $E_b/N_0=%g$ dB, $\\rho=%.1f$', ...
    ebn0, rho_sel), 'Interpreter','latex');

%% Diagnostic Plot 5: distance/step summary
figure; clf;
plot(Last2Summary_r.step, Last2Summary_r.P_tie, '-o'); hold on;
plot(Last2Summary_r.step, Last2Summary_r.P_full_swap_given_tie, '-s');
plot(Last2Summary_r.step, Last2Summary_r.P_tie_and_full_swap, '-^');
grid on;
xlabel('Bitonic step in last two stages','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
legend({'$P(\mathrm{MSB3\ tie})$', ...
        '$P(\mathrm{full\ swap}\mid\mathrm{MSB3\ tie})$', ...
        '$P(\mathrm{MSB3\ tie\ and\ full\ swap})$'}, ...
        'Interpreter','latex', 'Location','best');
title(sprintf('Tie and full-swap statistics by step, $E_b/N_0=%g$ dB, $\\rho=%.1f$', ...
    ebn0, rho_sel), 'Interpreter','latex');

%% Diagnostic Plot 6: full-swap fraction by last-stage distance
figure; clf;
bar(LastStageSummary_r.distance, LastStageSummary_r.fraction_of_last_stage_full_swaps);
grid on;
xlabel('Last-stage distance $j$','Interpreter','latex');
ylabel('Fraction of last-stage full swaps','Interpreter','latex');
title(sprintf('Distribution of MSB3-tie full-swaps by last-stage distance, $E_b/N_0=%g$ dB, $\\rho=%.1f$', ...
    ebn0, rho_sel), 'Interpreter','latex');

%% Diagnostic Plot 7: CAE concentration curve
figure; clf;
plot(1:height(T_last2_r), T_last2_r.cum_fraction, 'LineWidth', 1.5);
grid on;
xlabel('Top-$K$ CAEs ranked by full-swap count','Interpreter','latex');
ylabel('Cumulative fraction of full swaps','Interpreter','latex');
title(sprintf('CAE-level concentration of MSB3-tie full-swaps, $E_b/N_0=%g$ dB, $\\rho=%.1f$', ...
    ebn0, rho_sel), 'Interpreter','latex');
ylim([0 1]);