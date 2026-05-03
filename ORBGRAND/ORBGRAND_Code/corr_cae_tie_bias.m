% ==========================================================
% AR(1) correlated-noise sweep for MSB3 sorter swap behavior
%
% Measures:
%   1) P(full swap | MSB3 tie)
%   2) P(full keep | MSB3 tie)
%   3) P(actual MSB3 sorter swap) per comparison
%
% Requires updated MEX with 5 outputs:
%   [sorted_idx, tie_count, total_count, tie_swap_count, swap_count] = ...
%       corr_cae_tie_bias_mex(y_soft)
%
% Important:
%   swap_count = actual swaps performed by MSB3-only sorter
%   tie_swap_count = full-precision swaps only inside MSB3 ties
% ==========================================================

clear; clc;

%% User options
code_class = 'CRC';     % 'CAPOLAR' or 'CRC'
ebn0       = 9;
nFrames    = 3e4;

rho_list = [0 0.3 0.5 0.7 0.9 -0.3 -0.5 -0.7 -0.9];

rho_plot = 0.9;         % rho used for detailed diagnostic plots

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

P_tie_region                 = zeros(nrho,1);
P_full_swap_given_tie_region = zeros(nrho,1);
P_full_keep_given_tie_region = zeros(nrho,1);
P_tie_and_full_swap_region   = zeros(nrho,1);

P_actual_swap_region         = zeros(nrho,1);

P_full_swap_given_tie_step   = nan(n_steps,nrho);
P_tie_step                   = zeros(n_steps,nrho);
P_actual_swap_step           = zeros(n_steps,nrho);

%% Full per-CAE storage for diagnostics
tie_sum_all      = zeros(n_steps,n_cae,nrho);
total_sum_all    = zeros(n_steps,n_cae,nrho);
tie_swap_sum_all = zeros(n_steps,n_cae,nrho);
swap_sum_all     = zeros(n_steps,n_cae,nrho);

%% Main sweep
for rid = 1:nrho

    rho = rho_list(rid);

    fprintf('\n========================================\n');
    fprintf('Running rho = %.2f\n', rho);
    fprintf('========================================\n');

    tie_sum      = zeros(n_steps,n_cae);
    total_sum    = zeros(n_steps,n_cae);
    tie_swap_sum = zeros(n_steps,n_cae);
    swap_sum     = zeros(n_steps,n_cae);

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

        %% Correlation-aware BPSK LLR from Eq. (15)
        % Paper convention gives lambda = log P(x=-1|y)/P(x=+1|y),
        % while nrSymbolDemodulate gives positive LLR for +1.
        % Therefore we flip the sign to match the existing ORBGRAND convention.
        
        y = real(rSig(:)).';     % force real row vector
        y_prev = [y(1), y(1:end-1)];
        
        if abs(rho) < 1e-12
            y_soft = 2*y/sigma2;
        else
            lambda_paper = ( ...
                -2*y + 2*rho*y_prev ...
                + abs(y_prev - rho*y - rho) ...
                - abs(y_prev - rho*y + rho) ...
                ) ./ (sigma2*(1-rho^2));
        
            y_soft = -lambda_paper;
            y_soft(1) = 2*y(1)/sigma2;
        end

        y_soft = double(real(y_soft(:)));   % force real double column
        
        [~, tie_c, total_c, tie_swap_c, swap_c] = corr_cae_tie_bias_mex(y_soft);

        tie_sum      = tie_sum      + double(tie_c);
        total_sum    = total_sum    + double(total_c);
        tie_swap_sum = tie_swap_sum + double(tie_swap_c);
        swap_sum     = swap_sum     + double(swap_c);
    end

    %% Store full per-CAE counters
    tie_sum_all(:,:,rid)      = tie_sum;
    total_sum_all(:,:,rid)    = total_sum;
    tie_swap_sum_all(:,:,rid) = tie_swap_sum;
    swap_sum_all(:,:,rid)     = swap_sum;

    %% Region-level summary: last two stages
    total_ties_region         = sum(tie_sum(last2_steps,:),'all');
    total_full_swaps_region   = sum(tie_swap_sum(last2_steps,:),'all');
    total_actual_swaps_region = sum(swap_sum(last2_steps,:),'all');
    total_comp_region         = sum(total_sum(last2_steps,:),'all');

    P_tie_region(rid)                 = total_ties_region / total_comp_region;
    P_full_swap_given_tie_region(rid) = total_full_swaps_region / total_ties_region;
    P_full_keep_given_tie_region(rid) = 1 - P_full_swap_given_tie_region(rid);
    P_tie_and_full_swap_region(rid)   = total_full_swaps_region / total_comp_region;

    P_actual_swap_region(rid)         = total_actual_swaps_region / total_comp_region;

    %% Step-level summary
    for st = 1:n_steps
        total_ties_st         = sum(tie_sum(st,:),'all');
        total_full_swaps_st   = sum(tie_swap_sum(st,:),'all');
        total_actual_swaps_st = sum(swap_sum(st,:),'all');
        total_comp_st         = sum(total_sum(st,:),'all');

        P_tie_step(st,rid) = total_ties_st / total_comp_st;
        P_actual_swap_step(st,rid) = total_actual_swaps_st / total_comp_st;

        if total_ties_st > 0
            P_full_swap_given_tie_step(st,rid) = total_full_swaps_st / total_ties_st;
        end
    end

    fprintf('\nrho = %.2f summary, last two stages:\n', rho);
    fprintf('P(tie)                         = %.4f\n', P_tie_region(rid));
    fprintf('P(full swap | MSB tie)         = %.4f\n', P_full_swap_given_tie_region(rid));
    fprintf('P(full keep | MSB tie)         = %.4f\n', P_full_keep_given_tie_region(rid));
    fprintf('P(tie and full swap)           = %.4f\n', P_tie_and_full_swap_region(rid));
    fprintf('P(actual MSB3 sorter swap)     = %.4f\n', P_actual_swap_region(rid));
end

%% Summary table
Summary = table( ...
    rho_list(:), ...
    P_tie_region, ...
    P_full_swap_given_tie_region, ...
    P_full_keep_given_tie_region, ...
    P_tie_and_full_swap_region, ...
    P_actual_swap_region, ...
    'VariableNames', ...
    {'rho','P_tie','P_full_swap_given_MSB_tie', ...
     'P_full_keep_given_MSB_tie','P_tie_and_full_swap', ...
     'P_actual_MSB3_swap'} ...
);

disp(Summary);

save('AR1_MSB3_SWAP_SWEEP_WITH_DIAGNOSTICS.mat', ...
    'Summary', ...
    'P_full_swap_given_tie_step', 'P_tie_step', 'P_actual_swap_step', ...
    'tie_sum_all', 'total_sum_all', 'tie_swap_sum_all', 'swap_sum_all', ...
    'rho_list', 'step_outer', 'step_local', 'step_dist', ...
    'last2_steps', 'last_stage_steps', 'code_label', 'ebn0', ...
    'snr_db', 'sigma2');

%% ==========================================================
%% Sweep plots
%% ==========================================================

%% Plot 1: full swap/keep ratio inside MSB ties vs rho
figure; clf;
plot(rho_list, P_full_swap_given_tie_region, '-o', 'LineWidth', 1.5); hold on;
plot(rho_list, P_full_keep_given_tie_region, '-s', 'LineWidth', 1.5);
grid on;
xlabel('$\rho$','Interpreter','latex');
ylabel('Conditional probability','Interpreter','latex');
legend({'$P(\mathrm{full\ swap}\mid\mathrm{MSB3\ tie})$', ...
        '$P(\mathrm{full\ keep}\mid\mathrm{MSB3\ tie})$'}, ...
        'Interpreter','latex', 'Location','best');
title(sprintf('MSB3-tie resolution bias, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% Plot 2: tie probability and tie-and-full-swap vs rho
figure; clf;
plot(rho_list, P_tie_region, '-o', 'LineWidth', 1.5); hold on;
plot(rho_list, P_tie_and_full_swap_region, '-s', 'LineWidth', 1.5);
grid on;
xlabel('$\rho$','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
legend({'$P(\mathrm{MSB3\ tie})$', ...
        '$P(\mathrm{MSB3\ tie\ and\ full\ swap})$'}, ...
        'Interpreter','latex', 'Location','best');
title(sprintf('Tie occurrence under AR(1) correlated noise, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% Plot 3: actual MSB3 sorter swap probability vs rho
figure; clf;
plot(rho_list, P_actual_swap_region, '-o', 'LineWidth', 1.5);
grid on;
xlabel('$\rho$','Interpreter','latex');
ylabel('$P(\mathrm{actual\ MSB3\ swap})$','Interpreter','latex');
title(sprintf('Actual MSB3 sorter swap probability, last two stages, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% Plot 4: step-level full-swap-given-tie in last two stages
figure; clf;
hold on;

for rid = 1:nrho
    plot(last2_steps, P_full_swap_given_tie_step(last2_steps,rid), '-o', ...
        'DisplayName', sprintf('\\rho=%.1f', rho_list(rid)));
end

grid on;
xlabel('Bitonic step','Interpreter','latex');
ylabel('$P(\mathrm{full\ swap}\mid\mathrm{MSB3\ tie})$', ...
    'Interpreter','latex');
legend('Location','best');
title(sprintf('Step-level full-swap probability under MSB3 ties, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% Plot 5: actual swap probability per bitonic step
figure; clf;
hold on;

for rid = 1:nrho
    plot(1:n_steps, P_actual_swap_step(:,rid), '-o', ...
        'LineWidth', 1.3, ...
        'DisplayName', sprintf('\\rho=%.1f', rho_list(rid)));
end

grid on;
xlabel('Bitonic step','Interpreter','latex');
ylabel('$P(\mathrm{actual\ Full-precision\ swap})$','Interpreter','latex');
legend('Location','best');
title(sprintf('Actual Full-precision sorter swap probability per comparison, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% Plot 6: early-step actual swap probability
early_steps = 1:min(10,n_steps);

figure; clf;
hold on;

for rid = 1:nrho
    plot(early_steps, P_actual_swap_step(early_steps,rid), '-o', ...
        'LineWidth', 1.3, ...
        'DisplayName', sprintf('\\rho=%.1f', rho_list(rid)));
end

grid on;
xlabel('Early bitonic step','Interpreter','latex');
ylabel('$P(\mathrm{actual\ MSB3\ swap})$','Interpreter','latex');
legend('Location','best');
title(sprintf('Early-step actual swap probability under AR(1) noise, $E_b/N_0=%g$ dB', ebn0), ...
    'Interpreter','latex');

%% ==========================================================
%% Detailed diagnostics for selected rho
%% ==========================================================

[~, rid_plot] = min(abs(rho_list - rho_plot));
rho_sel = rho_list(rid_plot);

tie_sum_r      = tie_sum_all(:,:,rid_plot);
total_sum_r    = total_sum_all(:,:,rid_plot);
tie_swap_sum_r = tie_swap_sum_all(:,:,rid_plot);
swap_sum_r     = swap_sum_all(:,:,rid_plot);

p_tie_r = zeros(size(tie_sum_r));
valid_total = total_sum_r > 0;
p_tie_r(valid_total) = tie_sum_r(valid_total) ./ total_sum_r(valid_total);

p_full_swap_given_tie_r = nan(size(tie_sum_r));
valid_tie = tie_sum_r > 0;
p_full_swap_given_tie_r(valid_tie) = tie_swap_sum_r(valid_tie) ./ tie_sum_r(valid_tie);

p_tie_and_full_swap_r = zeros(size(tie_sum_r));
p_tie_and_full_swap_r(valid_total) = tie_swap_sum_r(valid_total) ./ total_sum_r(valid_total);

p_actual_swap_r = zeros(size(swap_sum_r));
p_actual_swap_r(valid_total) = swap_sum_r(valid_total) ./ total_sum_r(valid_total);

%% Last-two-stage step summary
rows = [];
total_last2_full_swaps   = sum(tie_swap_sum_r(last2_steps,:),'all');
total_last2_actual_swaps = sum(swap_sum_r(last2_steps,:),'all');

for st = last2_steps

    total_ties_st         = sum(tie_sum_r(st,:),'all');
    total_full_swaps_st   = sum(tie_swap_sum_r(st,:),'all');
    total_actual_swaps_st = sum(swap_sum_r(st,:),'all');
    total_comp_st         = sum(total_sum_r(st,:),'all');

    P_tie_st = total_ties_st / total_comp_st;
    P_full_swap_given_tie_st = total_full_swaps_st / total_ties_st;
    P_tie_and_full_swap_st = total_full_swaps_st / total_comp_st;
    P_actual_swap_st = total_actual_swaps_st / total_comp_st;

    frac_last2_full_swaps_st = total_full_swaps_st / total_last2_full_swaps;
    frac_last2_actual_swaps_st = total_actual_swaps_st / total_last2_actual_swaps;

    rows = [rows; ...
        st, step_outer(st), step_local(st), step_dist(st), ...
        total_ties_st, total_full_swaps_st, total_actual_swaps_st, total_comp_st, ...
        P_tie_st, P_full_swap_given_tie_st, P_tie_and_full_swap_st, P_actual_swap_st, ...
        frac_last2_full_swaps_st, frac_last2_actual_swaps_st];
end

Last2Summary_r = array2table(rows, ...
    'VariableNames', ...
    {'step','outer_stage','local_step','distance', ...
     'total_ties','full_swaps_inside_ties','actual_msb3_swaps','total_comparisons', ...
     'P_tie','P_full_swap_given_tie','P_tie_and_full_swap','P_actual_MSB3_swap', ...
     'fraction_of_last2_full_swaps','fraction_of_last2_actual_swaps'});

fprintf('\nLAST-TWO-STAGE SUMMARY FOR rho = %.2f:\n', rho_sel);
disp(Last2Summary_r);

%% Last-stage distance summary
rows = [];
total_last_stage_full_swaps   = sum(tie_swap_sum_r(last_stage_steps,:),'all');
total_last_stage_actual_swaps = sum(swap_sum_r(last_stage_steps,:),'all');

for st = last_stage_steps(:)'

    d = step_dist(st);

    total_ties_d         = sum(tie_sum_r(st,:),'all');
    total_full_swaps_d   = sum(tie_swap_sum_r(st,:),'all');
    total_actual_swaps_d = sum(swap_sum_r(st,:),'all');
    total_comp_d         = sum(total_sum_r(st,:),'all');

    P_tie_d = total_ties_d / total_comp_d;
    P_full_swap_given_tie_d = total_full_swaps_d / total_ties_d;
    P_tie_and_full_swap_d = total_full_swaps_d / total_comp_d;
    P_actual_swap_d = total_actual_swaps_d / total_comp_d;

    frac_last_stage_full_swaps_d = total_full_swaps_d / total_last_stage_full_swaps;
    frac_last_stage_actual_swaps_d = total_actual_swaps_d / total_last_stage_actual_swaps;

    rows = [rows; ...
        d, st, total_ties_d, total_full_swaps_d, total_actual_swaps_d, total_comp_d, ...
        P_tie_d, P_full_swap_given_tie_d, P_tie_and_full_swap_d, P_actual_swap_d, ...
        frac_last_stage_full_swaps_d, frac_last_stage_actual_swaps_d];
end

LastStageSummary_r = array2table(rows, ...
    'VariableNames', ...
    {'distance','step','total_ties','full_swaps_inside_ties','actual_msb3_swaps','total_comparisons', ...
     'P_tie','P_full_swap_given_tie','P_tie_and_full_swap','P_actual_MSB3_swap', ...
     'fraction_of_last_stage_full_swaps','fraction_of_last_stage_actual_swaps'});

LastStageSummary_r = sortrows(LastStageSummary_r, 'distance', 'ascend');

fprintf('\nLAST-STAGE DISTANCE SUMMARY FOR rho = %.2f:\n', rho_sel);
disp(LastStageSummary_r);

%% Uniformity check
valid_last2_psgt = p_full_swap_given_tie_r(last2_steps,:);
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

%% Correlation checks
corr_full_last2 = corr(Last2Summary_r.P_tie, Last2Summary_r.fraction_of_last2_full_swaps);
corr_full_last_stage = corr(LastStageSummary_r.P_tie, LastStageSummary_r.fraction_of_last_stage_full_swaps);

corr_actual_last2 = corr(Last2Summary_r.P_actual_MSB3_swap, Last2Summary_r.fraction_of_last2_actual_swaps);
corr_actual_last_stage = corr(LastStageSummary_r.P_actual_MSB3_swap, LastStageSummary_r.fraction_of_last_stage_actual_swaps);

fprintf('\nDISTANCE EFFECT CHECK, rho = %.2f:\n', rho_sel);
fprintf('corr(P_tie, full-swap fraction) last two stages       = %.4f\n', corr_full_last2);
fprintf('corr(P_tie, full-swap fraction) last stage            = %.4f\n', corr_full_last_stage);
fprintf('corr(P_actual_swap, actual-swap fraction) last two    = %.4f\n', corr_actual_last2);
fprintf('corr(P_actual_swap, actual-swap fraction) last stage  = %.4f\n', corr_actual_last_stage);

%% Diagnostic Plot 7: heatmap full swap inside ties
figure; clf;
imagesc(p_full_swap_given_tie_r(last2_steps,:));
colorbar;
caxis([0 1]);
xlabel('CAE index','Interpreter','latex');
ylabel('Last-two-stage step','Interpreter','latex');
yticks(1:numel(last2_steps));
yticklabels(last2_dist_labels);
title(sprintf('$P(\\mathrm{full\\ swap}\\mid\\mathrm{MSB3\\ tie})$, $E_b/N_0=%g$ dB, $\\rho=%.1f$', ...
    ebn0, rho_sel), 'Interpreter','latex');

%% Diagnostic Plot 8: heatmap actual MSB3 swap probability
figure; clf;
imagesc(p_actual_swap_r(last2_steps,:));
colorbar;
caxis([0 1]);
xlabel('CAE index','Interpreter','latex');
ylabel('Last-two-stage step','Interpreter','latex');
yticks(1:numel(last2_steps));
yticklabels(last2_dist_labels);
title(sprintf('$P(\\mathrm{actual\\ MSB3\\ swap})$, $E_b/N_0=%g$ dB, $\\rho=%.1f$', ...
    ebn0, rho_sel), 'Interpreter','latex');

%% Diagnostic Plot 9: step summary
figure; clf;
plot(Last2Summary_r.step, Last2Summary_r.P_tie, '-o'); hold on;
plot(Last2Summary_r.step, Last2Summary_r.P_full_swap_given_tie, '-s');
plot(Last2Summary_r.step, Last2Summary_r.P_tie_and_full_swap, '-^');
plot(Last2Summary_r.step, Last2Summary_r.P_actual_MSB3_swap, '-d');
grid on;
xlabel('Bitonic step in last two stages','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
legend({'$P(\mathrm{MSB3\ tie})$', ...
        '$P(\mathrm{full\ swap}\mid\mathrm{MSB3\ tie})$', ...
        '$P(\mathrm{MSB3\ tie\ and\ full\ swap})$', ...
        '$P(\mathrm{actual\ MSB3\ swap})$'}, ...
        'Interpreter','latex', 'Location','best');
title(sprintf('Tie and swap statistics by step, $E_b/N_0=%g$ dB, $\\rho=%.1f$', ...
    ebn0, rho_sel), 'Interpreter','latex');

%% Diagnostic Plot 10: actual swap fraction by last-stage distance
figure; clf;
bar(LastStageSummary_r.distance, LastStageSummary_r.fraction_of_last_stage_actual_swaps);
grid on;
xlabel('Last-stage distance $j$','Interpreter','latex');
ylabel('Fraction of last-stage actual MSB3 swaps','Interpreter','latex');
title(sprintf('Distribution of actual MSB3 swaps by last-stage distance, $E_b/N_0=%g$ dB, $\\rho=%.1f$', ...
    ebn0, rho_sel), 'Interpreter','latex');