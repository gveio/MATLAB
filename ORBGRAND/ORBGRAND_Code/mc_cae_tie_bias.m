% ==========================================================
% PA lost-swap diagnostics for MSB3 ties
%
% Supports these conclusions:
%   1) No CAE-level structure
%   2) No swap-bias to exploit
%   3) Only distance/tie probability matters
%   4) PA tie-breaking is stochastic/adaptive
% ==========================================================

clear; clc;

%% User options
code_class = 'CRC';     % 'CAPOLAR' or 'CRC'
ebn0       = 7;
nFrames    = 3e4;

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

%% Step metadata
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

%% Regions
last_stage_steps = find(step_outer == m_bits);

last2_outer_start = (m_bits-2)*(m_bits-1)/2 + 1;
last2_steps       = last2_outer_start:n_steps;

last2_dist_labels = strings(numel(last2_steps),1);
for i = 1:numel(last2_steps)
    st = last2_steps(i);
    last2_dist_labels(i) = sprintf('step %d, d=%d', st, step_dist(st));
end

fprintf('n = %d, n_effective = %d\n', n, n_effective);
fprintf('Total bitonic steps = %d\n', n_steps);
fprintf('Last stage steps: ');
fprintf('%d ', last_stage_steps);
fprintf('\n');
fprintf('Last two outer stages: steps %d to %d\n', last2_steps(1), last2_steps(end));

%% Accumulators
tie_sum      = zeros(n_steps, n_cae, numSNR);
total_sum    = zeros(n_steps, n_cae, numSNR);
tie_swap_sum = zeros(n_steps, n_cae, numSNR);

%% Monte Carlo
for ii = 1:numSNR

    sigma2 = 1/(10^(0.1*snr_db(ii)));
    awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);

    fprintf('\nRunning %s Eb/N0=%.2f dB\n', code_class, ebn0(ii));

    for frm = 1:nFrames

        if mod(frm,1000) == 0
            fprintf('  Frame %d / %d\n', frm, nFrames);
        end

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

%% ==========================================================
%% Compute probabilities
%% ==========================================================

p_tie = zeros(size(tie_sum));
valid_total = total_sum > 0;
p_tie(valid_total) = tie_sum(valid_total) ./ total_sum(valid_total);

p_swap_given_tie = nan(size(tie_sum));
valid_tie = tie_sum > 0;
p_swap_given_tie(valid_tie) = tie_swap_sum(valid_tie) ./ tie_sum(valid_tie);

p_tie_and_swap = zeros(size(tie_sum));
p_tie_and_swap(valid_total) = tie_swap_sum(valid_total) ./ total_sum(valid_total);

%% ==========================================================
%% Compact terminal diagnostics
%% ==========================================================

ii = 1;

fprintf('\n==========================================================\n');
fprintf('KEY DIAGNOSTICS, Eb/N0 = %.2f dB\n', ebn0(ii));
fprintf('==========================================================\n');

%% Region summary
regions = {
    'LAST_STAGE',     last_stage_steps;
    'LAST_TWO_STAGE', last2_steps
};

region_rows = [];

for rr = 1:size(regions,1)

    region_steps = regions{rr,2};

    total_ties_region  = sum(tie_sum(region_steps,:,ii),'all');
    total_swaps_region = sum(tie_swap_sum(region_steps,:,ii),'all');
    total_comp_region  = sum(total_sum(region_steps,:,ii),'all');

    P_tie_region = total_ties_region / total_comp_region;
    P_swap_given_tie_region = total_swaps_region / total_ties_region;
    P_tie_and_swap_region = total_swaps_region / total_comp_region;

    region_rows = [region_rows; ...
        rr, total_ties_region, total_swaps_region, total_comp_region, ...
        P_tie_region, P_swap_given_tie_region, P_tie_and_swap_region];
end

RegionSummary = array2table(region_rows, ...
    'VariableNames', ...
    {'region_id','total_ties','PA_swaps','total_comparisons', ...
     'P_tie','P_PA_swap_given_tie','P_tie_and_PA_swap'});

RegionSummary.region = ["LAST_STAGE"; "LAST_TWO_STAGE"];
RegionSummary = movevars(RegionSummary, 'region', 'After', 'region_id');

fprintf('\nREGION SUMMARY:\n');
disp(RegionSummary);

%% Last-two-stage step/distance summary
rows = [];
total_last2_swaps = sum(tie_swap_sum(last2_steps,:,ii),'all');

for st = last2_steps

    total_ties_st  = sum(tie_sum(st,:,ii),'all');
    total_swaps_st = sum(tie_swap_sum(st,:,ii),'all');
    total_comp_st  = sum(total_sum(st,:,ii),'all');

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

Last2Summary = array2table(rows, ...
    'VariableNames', ...
    {'step','outer_stage','local_step','distance', ...
     'total_ties','PA_swaps','total_comparisons', ...
     'P_tie','P_PA_swap_given_tie','P_tie_and_PA_swap', ...
     'fraction_of_last2_PA_swaps'});

fprintf('\nLAST-TWO-STAGE STEP/DISTANCE SUMMARY:\n');
disp(Last2Summary);

%% Last-stage distance summary
rows = [];
total_last_stage_swaps = sum(tie_swap_sum(last_stage_steps,:,ii),'all');

for st = last_stage_steps(:)'

    d = step_dist(st);

    total_ties_d  = sum(tie_sum(st,:,ii),'all');
    total_swaps_d = sum(tie_swap_sum(st,:,ii),'all');
    total_comp_d  = sum(total_sum(st,:,ii),'all');

    P_tie_d = total_ties_d / total_comp_d;
    P_swap_given_tie_d = total_swaps_d / total_ties_d;
    P_tie_and_swap_d = total_swaps_d / total_comp_d;
    frac_last_stage_swaps_d = total_swaps_d / total_last_stage_swaps;

    rows = [rows; ...
        d, st, total_ties_d, total_swaps_d, total_comp_d, ...
        P_tie_d, P_swap_given_tie_d, P_tie_and_swap_d, ...
        frac_last_stage_swaps_d];
end

LastStageSummary = array2table(rows, ...
    'VariableNames', ...
    {'distance','step','total_ties','PA_swaps','total_comparisons', ...
     'P_tie','P_PA_swap_given_tie','P_tie_and_PA_swap', ...
     'fraction_of_last_stage_PA_swaps'});

LastStageSummary = sortrows(LastStageSummary, 'distance', 'ascend');

fprintf('\nLAST-STAGE DISTANCE SUMMARY:\n');
disp(LastStageSummary);

fprintf('\nDistance-wise PA swap contribution in LAST STAGE:\n');
for krow = 1:height(LastStageSummary)
    fprintf('d=%3d : fraction=%.4f, P(tie)=%.4f, P(PA swap|tie)=%.4f, P(tie and PA swap)=%.4e\n', ...
        LastStageSummary.distance(krow), ...
        LastStageSummary.fraction_of_last_stage_PA_swaps(krow), ...
        LastStageSummary.P_tie(krow), ...
        LastStageSummary.P_PA_swap_given_tie(krow), ...
        LastStageSummary.P_tie_and_PA_swap(krow));
end

%% Uniformity of P(PA swap | tie)
valid_last2_psgt = p_swap_given_tie(last2_steps,:,ii);
valid_last2_psgt = valid_last2_psgt(~isnan(valid_last2_psgt));

mean_psgt = mean(valid_last2_psgt);
std_psgt  = std(valid_last2_psgt);
min_psgt  = min(valid_last2_psgt);
max_psgt  = max(valid_last2_psgt);

fprintf('\nUNIFORMITY OF P(PA swap | tie) IN LAST TWO STAGES:\n');
fprintf('mean = %.4f\n', mean_psgt);
fprintf('std  = %.4f\n', std_psgt);
fprintf('min  = %.4f\n', min_psgt);
fprintf('max  = %.4f\n', max_psgt);

%% Correlation: PA swap fraction follows tie probability
corr_last2 = corr(Last2Summary.P_tie, Last2Summary.fraction_of_last2_PA_swaps);
corr_last_stage = corr(LastStageSummary.P_tie, LastStageSummary.fraction_of_last_stage_PA_swaps);

fprintf('\nDISTANCE EFFECT CHECK:\n');
fprintf('corr(P_tie, PA-swap fraction) last two stages = %.4f\n', corr_last2);
fprintf('corr(P_tie, PA-swap fraction) last stage      = %.4f\n', corr_last_stage);

%% CAE concentration ranking: last two stages
last2_mask = false(n_steps,n_cae);
last2_mask(last2_steps,:) = true;

valid = last2_mask & (tie_swap_sum(:,:,ii) > 0);

[step_id, cae_id] = find(valid);

T_last2 = table( ...
    step_id, ...
    step_outer(step_id), ...
    step_local(step_id), ...
    step_dist(step_id), ...
    cae_id, ...
    tie_sum(valid), ...
    tie_swap_sum(valid), ...
    total_sum(valid), ...
    p_tie(valid), ...
    p_swap_given_tie(valid), ...
    p_tie_and_swap(valid), ...
    'VariableNames', ...
    {'step','outer_stage','local_step','distance','cae', ...
     'tie_count','PA_swap_count','total_count', ...
     'P_tie','P_PA_swap_given_tie','P_tie_and_PA_swap'} ...
);

T_last2 = sortrows(T_last2, 'PA_swap_count', 'descend');

total_swaps_last2 = sum(tie_swap_sum(last2_steps,:,ii),'all');

T_last2.cum_PA_swaps = cumsum(T_last2.PA_swap_count);
T_last2.cum_fraction = T_last2.cum_PA_swaps / total_swaps_last2;

fprintf('\nCAE CONCENTRATION IN LAST TWO STAGES:\n');
fprintf('Total CAEs in last two stages : %d\n', numel(last2_steps)*n_cae);
fprintf('CAEs with PA swaps            : %d\n', height(T_last2));
fprintf('Total PA swaps                : %d\n', total_swaps_last2);

Klist = [10 50 100 200 400 800 1200 1600];

fprintf('\nTop-K cumulative PA swap fraction:\n');
for kk = 1:length(Klist)
    Ksel = Klist(kk);
    if Ksel <= height(T_last2)
        fprintf('Top %4d -> fraction=%.4f\n', Ksel, T_last2.cum_fraction(Ksel));
    end
end

%% Terminal conclusion
fprintf('\n==========================================================\n');
fprintf('INTERPRETATION CHECKLIST\n');
fprintf('==========================================================\n');

fprintf('1) No CAE-level structure:\n');
fprintf('   If Top-K curve rises slowly, swaps are distributed across CAEs.\n');

fprintf('2) No bias to exploit:\n');
fprintf('   P(PA swap | tie) is nearly constant; mean %.4f, std %.4f.\n', ...
    mean_psgt, std_psgt);

fprintf('3) Only distance/tie probability matters:\n');
fprintf('   PA-swap fraction correlates with P(tie): last2 corr %.4f.\n', corr_last2);

fprintf('4) PA is stochastic/adaptive:\n');
fprintf('   Since P(PA swap | tie) ≈ %.2f, fixed-swap is wrong most tie cases.\n', mean_psgt);
fprintf('   Therefore fixed CAE masks are weak; distance-limited PA is the sensible hardware experiment.\n');

%% ==========================================================
%% Store
%% ==========================================================

result = struct();

result.p_tie             = p_tie;
result.p_swap_given_tie  = p_swap_given_tie;
result.p_tie_and_swap    = p_tie_and_swap;

result.tie_sum           = tie_sum;
result.total_sum         = total_sum;
result.tie_swap_sum      = tie_swap_sum;

result.step_outer        = step_outer;
result.step_local        = step_local;
result.step_dist         = step_dist;

result.last_stage_steps  = last_stage_steps;
result.last2_steps       = last2_steps;

result.RegionSummary     = RegionSummary;
result.Last2Summary      = Last2Summary;
result.LastStageSummary  = LastStageSummary;
result.T_last2_ranked    = T_last2;

result.mean_psgt_last2   = mean_psgt;
result.std_psgt_last2    = std_psgt;
result.corr_last2        = corr_last2;
result.corr_last_stage   = corr_last_stage;

result.ebn0              = ebn0;
result.n_steps           = n_steps;
result.n_cae             = n_cae;
result.code_label        = code_label;

save('PA_LOST_SWAP_DIAGNOSTICS_MINIMAL.mat','result');

%% ==========================================================
%% Minimum plots
%% ==========================================================

%% Plot 1: P(PA swap | tie) heatmap
% Shows: no clear CAE-level/bias structure; values are nearly uniform.
figure; clf;
imagesc(p_swap_given_tie(last2_steps,:,ii));
colorbar;
caxis([0 1]);
xlabel('CAE index','Interpreter','latex');
ylabel('Last-two-stage step','Interpreter','latex');
yticks(1:numel(last2_steps));
yticklabels(last2_dist_labels);
title(sprintf('$P(\\mathrm{PA\\ swap}\\mid\\mathrm{MSB\\ tie})$, $E_b/N_0=%g$ dB', ebn0(ii)), ...
    'Interpreter','latex');

%% Plot 2: distance/step summary
% Shows: P(PA swap | tie) nearly flat, P(tie) changes strongly with distance.
figure; clf;
plot(Last2Summary.step, Last2Summary.P_tie, '-o'); hold on;
plot(Last2Summary.step, Last2Summary.P_PA_swap_given_tie, '-s');
plot(Last2Summary.step, Last2Summary.P_tie_and_PA_swap, '-^');
grid on;
xlabel('Bitonic step in last two stages','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
legend({'$P(\mathrm{tie})$', ...
        '$P(\mathrm{PA\ swap}\mid\mathrm{tie})$', ...
        '$P(\mathrm{tie\ and\ PA\ swap})$'}, ...
        'Interpreter','latex', 'Location','best');
title(sprintf('Tie and PA-swap statistics by step, $E_b/N_0=%g$ dB', ebn0(ii)), ...
    'Interpreter','latex');

%% Plot 3: PA swap fraction by last-stage distance
% Shows: only distance/tie probability matters.
figure; clf;
bar(LastStageSummary.distance, LastStageSummary.fraction_of_last_stage_PA_swaps);
grid on;
xlabel('Last-stage distance $j$','Interpreter','latex');
ylabel('Fraction of last-stage PA swaps','Interpreter','latex');
title(sprintf('Distribution of PA tie-swaps by last-stage distance, $E_b/N_0=%g$ dB', ebn0(ii)), ...
    'Interpreter','latex');

%% Plot 4: CAE concentration curve
% Shows: no small subset of CAEs dominates.
figure; clf;
plot(1:height(T_last2), T_last2.cum_fraction, 'LineWidth', 1.5);
grid on;
xlabel('Top-$K$ CAEs ranked by PA swap count','Interpreter','latex');
ylabel('Cumulative fraction of PA swaps','Interpreter','latex');
title(sprintf('CAE-level concentration of PA swaps, $E_b/N_0=%g$ dB', ebn0(ii)), ...
    'Interpreter','latex');
ylim([0 1]);