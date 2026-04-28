% ==========================================================
% PA2 useful-swap diagnostics for MSB3 ties
%
% Requires:
%   pa2_useful_swap_trace_mex.c
%
% MEX call:
% [idx_msb3, idx_pa2, tie_c, total_c, tie_swap_c, ...
%  swap_step, swap_cae, swap_a, swap_b] = pa2_useful_swap_trace_mex(y_soft)
%
% Goal:
%   Locate PA2 swaps in the final two stages that affect useful
%   final ranking positions: top-16, top-32, top-64, top-LW_MAX.
% ==========================================================

clear; clc;

%% User options
code_class = 'CRC';     % 'CAPOLAR' or 'CRC' or 'LDPC'
ebn0       = 7;
nFrames    = 3e4;

modulation = 'BPSK';
nmodbits   = 1;

LW_MAX = 104;
Qlist = [16 32 64 LW_MAX];

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

    case 'LDPC'
        S = open('ldpc_1024_512.mat');
        G = S.G; H = S.H;
        [k,n] = size(G);
        code_label = 'LDPC';

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
fprintf('Last two outer stages: steps %d to %d\n', last2_steps(1), last2_steps(end));

%% Accumulators
tie_sum      = zeros(n_steps, n_cae, numSNR);
total_sum    = zeros(n_steps, n_cae, numSNR);
tie_swap_sum = zeros(n_steps, n_cae, numSNR);

% Useful swap accumulators
nQ = length(Qlist);

useful_any_sum       = zeros(n_steps,n_cae,nQ,numSNR);
boundary_cross_sum   = zeros(n_steps,n_cae,nQ,numSNR);
relative_change_sum  = zeros(n_steps,n_cae,nQ,numSNR);
top_involved_sum     = zeros(n_steps,n_cae,nQ,numSNR);

total_traced_swaps = zeros(numSNR,1);
topQ_diff_frames   = zeros(nQ,numSNR);

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

        [idx_msb3, idx_pa2, tie_c, total_c, tie_swap_c, ...
         swap_step, swap_cae, swap_a, swap_b] = pa2_useful_swap_trace_mex(y_soft);

        tie_sum(:,:,ii)      = tie_sum(:,:,ii)      + double(tie_c);
        total_sum(:,:,ii)    = total_sum(:,:,ii)    + double(total_c);
        tie_swap_sum(:,:,ii) = tie_swap_sum(:,:,ii) + double(tie_swap_c);

        total_traced_swaps(ii) = total_traced_swaps(ii) + length(swap_step);

        %% Build rank maps
        rank_msb3 = zeros(1,n);
        rank_pa2  = zeros(1,n);

        for r = 1:n
            rank_msb3(double(idx_msb3(r))) = r;
            rank_pa2(double(idx_pa2(r)))   = r;
        end

        %% Frame-level top-Q difference
        for qid = 1:nQ
            Q = min(Qlist(qid), n);

            top_msb3 = sort(double(idx_msb3(1:Q)));
            top_pa2  = sort(double(idx_pa2(1:Q)));

            if ~isequal(top_msb3, top_pa2)
                topQ_diff_frames(qid,ii) = topQ_diff_frames(qid,ii) + 1;
            end
        end

        %% Classify traced PA2 swaps
        for e = 1:length(swap_step)

            st  = double(swap_step(e));
            cae = double(swap_cae(e));

            a = double(swap_a(e));
            b = double(swap_b(e));

            for qid = 1:nQ

                Q = min(Qlist(qid), n);

                a_top_msb3 = rank_msb3(a) <= Q;
                b_top_msb3 = rank_msb3(b) <= Q;

                a_top_pa2  = rank_pa2(a) <= Q;
                b_top_pa2  = rank_pa2(b) <= Q;

                top_involved = ...
                    a_top_msb3 || b_top_msb3 || ...
                    a_top_pa2  || b_top_pa2;

                boundary_cross = ...
                    xor(a_top_msb3, a_top_pa2) || ...
                    xor(b_top_msb3, b_top_pa2);

                relative_change = ...
                    sign(rank_msb3(a) - rank_msb3(b)) ~= ...
                    sign(rank_pa2(a)  - rank_pa2(b));

                useful_any = top_involved && ...
                             (boundary_cross || relative_change);

                if top_involved
                    top_involved_sum(st,cae,qid,ii) = top_involved_sum(st,cae,qid,ii) + 1;
                end

                if boundary_cross
                    boundary_cross_sum(st,cae,qid,ii) = boundary_cross_sum(st,cae,qid,ii) + 1;
                end

                if relative_change && top_involved
                    relative_change_sum(st,cae,qid,ii) = relative_change_sum(st,cae,qid,ii) + 1;
                end

                if useful_any
                    useful_any_sum(st,cae,qid,ii) = useful_any_sum(st,cae,qid,ii) + 1;
                end
            end
        end
    end
end

%% ==========================================================
%% Probability metrics
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
%% Terminal diagnostics
%% ==========================================================

ii = 1;

fprintf('\n==========================================================\n');
fprintf('PA2 USEFUL-SWAP DIAGNOSTICS, Eb/N0 = %.2f dB\n', ebn0(ii));
fprintf('==========================================================\n');

fprintf('Total traced PA2 swaps in last two stages: %d\n', total_traced_swaps(ii));

fprintf('\nFrame-level top-Q set difference between MSB3 and PA2:\n');
for qid = 1:nQ
    fprintf('Top-%d differs in %d / %d frames (%.4f)\n', ...
        Qlist(qid), ...
        topQ_diff_frames(qid,ii), ...
        nFrames, ...
        topQ_diff_frames(qid,ii)/nFrames);
end

%% Last-two-stage summary
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
     'total_ties','PA2_swaps','total_comparisons', ...
     'P_tie','P_PA2_swap_given_tie','P_tie_and_PA2_swap', ...
     'fraction_of_last2_PA2_swaps'});

fprintf('\nLAST-TWO-STAGE PA2 SWAP SUMMARY:\n');
disp(Last2Summary);

%% Useful-swap summary by Q
for qid = 1:nQ

    Q = Qlist(qid);

    useful_last2  = useful_any_sum(last2_steps,:,qid,ii);
    cross_last2   = boundary_cross_sum(last2_steps,:,qid,ii);
    rel_last2     = relative_change_sum(last2_steps,:,qid,ii);
    topinv_last2  = top_involved_sum(last2_steps,:,qid,ii);

    fprintf('\n==========================================================\n');
    fprintf('USEFUL PA2 SWAPS AFFECTING TOP-%d\n', Q);
    fprintf('==========================================================\n');

    fprintf('Top-involved PA2 swaps      : %d\n', sum(topinv_last2,'all'));
    fprintf('Boundary-cross PA2 swaps    : %d\n', sum(cross_last2,'all'));
    fprintf('Relative-change PA2 swaps   : %d\n', sum(rel_last2,'all'));
    fprintf('Useful-any PA2 swaps        : %d\n', sum(useful_last2,'all'));

    fprintf('Fraction useful among traced PA2 swaps: %.6f\n', ...
        sum(useful_last2,'all') / total_traced_swaps(ii));

    %% Rank useful CAEs
    useful_full = useful_any_sum(:,:,qid,ii);

    valid = false(n_steps,n_cae);
    valid(last2_steps,:) = useful_full(last2_steps,:) > 0;

    [step_id, cae_id] = find(valid);

    T_useful = table( ...
        step_id, ...
        step_outer(step_id), ...
        step_local(step_id), ...
        step_dist(step_id), ...
        cae_id, ...
        useful_full(valid), ...
        tie_swap_sum(valid), ...
        p_swap_given_tie(valid), ...
        p_tie_and_swap(valid), ...
        'VariableNames', ...
        {'step','outer_stage','local_step','distance','cae', ...
         'useful_count','PA2_swap_count','P_PA2_swap_given_tie','P_tie_and_PA2_swap'} ...
    );

    T_useful = sortrows(T_useful, 'useful_count', 'descend');

    total_useful_Q = sum(useful_last2,'all');

    if total_useful_Q > 0
        T_useful.cum_useful = cumsum(T_useful.useful_count);
        T_useful.cum_fraction = T_useful.cum_useful / total_useful_Q;
    end

    fprintf('\nTop useful CAEs for top-%d:\n', Q);
    disp(T_useful(1:min(20,height(T_useful)),:));

    if total_useful_Q > 0
        Klist = [10 50 100 200 400 800 1200 1600];

        fprintf('\nTop-K cumulative useful PA2 swap fraction, top-%d:\n', Q);
        for kk = 1:length(Klist)
            Ksel = Klist(kk);
            if Ksel <= height(T_useful)
                fprintf('Top %4d -> fraction=%.4f\n', ...
                    Ksel, T_useful.cum_fraction(Ksel));
            end
        end
    end
end

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

result.useful_any_sum      = useful_any_sum;
result.boundary_cross_sum  = boundary_cross_sum;
result.relative_change_sum = relative_change_sum;
result.top_involved_sum    = top_involved_sum;

result.total_traced_swaps = total_traced_swaps;
result.topQ_diff_frames   = topQ_diff_frames;
result.Qlist              = Qlist;

result.step_outer        = step_outer;
result.step_local        = step_local;
result.step_dist         = step_dist;

result.last_stage_steps  = last_stage_steps;
result.last2_steps       = last2_steps;

result.Last2Summary      = Last2Summary;

result.ebn0              = ebn0;
result.n_steps           = n_steps;
result.n_cae             = n_cae;
result.code_label        = code_label;

save('PA2_USEFUL_SWAP_DIAGNOSTICS.mat','result');

%% ==========================================================
%% Minimum plots
%% ==========================================================

ii = 1;

%% Plot 1: PA2 swap probability given MSB tie
figure; clf;
imagesc(p_swap_given_tie(last2_steps,:,ii));
colorbar;
caxis([0 1]);
xlabel('CAE index','Interpreter','latex');
ylabel('Last-two-stage step','Interpreter','latex');
yticks(1:numel(last2_steps));
yticklabels(last2_dist_labels);
title(sprintf('$P(\\mathrm{PA2\\ swap}\\mid\\mathrm{MSB\\ tie})$, $E_b/N_0=%g$ dB', ebn0(ii)), ...
    'Interpreter','latex');

%% Plot 2: Useful PA2 swaps affecting top-LWMAX
[~, lwid] = min(abs(Qlist - LW_MAX));

figure; clf;
imagesc(useful_any_sum(last2_steps,:,lwid,ii));
colorbar;
xlabel('CAE index','Interpreter','latex');
ylabel('Last-two-stage step','Interpreter','latex');
yticks(1:numel(last2_steps));
yticklabels(last2_dist_labels);
title(sprintf('Useful PA2 swaps affecting top-%d, $E_b/N_0=%g$ dB', Qlist(lwid), ebn0(ii)), ...
    'Interpreter','latex');

%% Plot 3: Step-level useful swaps for each Q
figure; clf;
hold on;

for qid = 1:nQ
    y = squeeze(sum(useful_any_sum(last2_steps,:,qid,ii),2));
    plot(last2_steps, y, '-o');
end

grid on;
xlabel('Bitonic step','Interpreter','latex');
ylabel('Useful PA2 swap count','Interpreter','latex');
legend(arrayfun(@(q) sprintf('top-%d',q), Qlist, 'UniformOutput', false), ...
    'Location','best');
title(sprintf('Useful PA2 swaps by step, $E_b/N_0=%g$ dB', ebn0(ii)), ...
    'Interpreter','latex');

%% Plot 4: CAE concentration of useful swaps for top-LWMAX
useful_full = useful_any_sum(:,:,lwid,ii);
valid = false(n_steps,n_cae);
valid(last2_steps,:) = useful_full(last2_steps,:) > 0;

[step_id, cae_id] = find(valid);

Tplot = table( ...
    step_id, ...
    cae_id, ...
    useful_full(valid), ...
    'VariableNames', {'step','cae','useful_count'} ...
);

Tplot = sortrows(Tplot, 'useful_count', 'descend');

if ~isempty(Tplot)
    Tplot.cum_useful = cumsum(Tplot.useful_count);
    Tplot.cum_fraction = Tplot.cum_useful / sum(Tplot.useful_count);

    figure; clf;
    plot(1:height(Tplot), Tplot.cum_fraction, 'LineWidth', 1.5);
    grid on;
    xlabel('Top-$K$ CAEs ranked by useful PA2 swap count','Interpreter','latex');
    ylabel('Cumulative fraction of useful PA2 swaps','Interpreter','latex');
    title(sprintf('Concentration of useful PA2 swaps, top-%d, $E_b/N_0=%g$ dB', Qlist(lwid), ebn0(ii)), ...
        'Interpreter','latex');
    ylim([0 1]);
end