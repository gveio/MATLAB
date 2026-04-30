% ==========================================================
% AR(1) correlation sweep for MSB3-tie full-precision behavior
%
% Measures:
%   P(full swap | MSB3 tie)
%   P(full keep | MSB3 tie)
%
% Goal:
%   Check if AWGN 40/60 swap/keep changes under correlated LLRs.
% ==========================================================

clear; clc;

%% User options
code_class = 'CRC';     % 'CAPOLAR' or 'CRC'
ebn0       = 7;
nFrames    = 3e4;

rho_list = [0 0.3 0.5 0.7 0.9 -0.3 -0.5 -0.7 -0.9];

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
awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);

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

fprintf('n = %d, n_effective = %d, steps = %d\n', n, n_effective, n_steps);
fprintf('Last two stages: steps %d to %d\n', last2_steps(1), last2_steps(end));

%% Results
nrho = length(rho_list);

P_tie_region            = zeros(nrho,1);
P_swap_given_tie_region = zeros(nrho,1);
P_keep_given_tie_region = zeros(nrho,1);
P_tie_and_swap_region   = zeros(nrho,1);

P_swap_step = nan(n_steps,nrho);
P_tie_step  = zeros(n_steps,nrho);

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

        %% Standard AWGN channel
        modOut = nrSymbolModulate([c mod_pad]', modulation);
        rSig   = awgnchan(modOut);

        y_soft_iid = nrSymbolDemodulate(rSig, modulation, sigma2);
        y_soft_iid = y_soft_iid(1:n)';

        %% AR(1)-correlated LLR sequence
        if rho == 0
            y_soft = y_soft_iid;
        else
            y_soft = zeros(size(y_soft_iid));
            y_soft(1) = y_soft_iid(1);

            for jj = 2:n
                y_soft(jj) = rho*y_soft(jj-1) + sqrt(1-rho^2)*y_soft_iid(jj);
            end
        end

        %% MEX diagnostic
        [~, tie_c, total_c, tie_swap_c] = pa_cae_tie_bias_mex(y_soft);

        tie_sum      = tie_sum      + double(tie_c);
        total_sum    = total_sum    + double(total_c);
        tie_swap_sum = tie_swap_sum + double(tie_swap_c);
    end

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

%% Table
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

save('AR1_MSB3_TIE_SWAP_SWEEP.mat', ...
    'Summary', 'P_swap_step', 'P_tie_step', ...
    'rho_list', 'step_outer', 'step_local', 'step_dist', ...
    'last2_steps', 'last_stage_steps', 'code_label', 'ebn0');

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
title(sprintf('Tie occurrence under AR(1) LLR correlation, $E_b/N_0=%g$ dB', ebn0), ...
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