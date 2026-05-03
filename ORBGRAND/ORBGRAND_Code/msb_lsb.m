clear; clc; close all;

% --- LaTeX style ---
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% Parameters
EbN0_list = [4 5 6 7];
n = 256;
k = 240;
R = k/n;

B_mag = 5;
LLR_max = 31;
N = 2e6;

bit_labels = {'MSB$_4$','MSB$_3$','MSB$_2$','LSB$_1$','LSB$_0$'};
P01_all = zeros(length(EbN0_list), B_mag);

random_ref = 0.25;
tol = 0.015;

%% Monte Carlo
for e = 1:length(EbN0_list)

    EbN0_dB = EbN0_list(e);

    snr_db = EbN0_dB + 10*log10(R);
    sigma2 = 1/(10^(0.1*snr_db));
    sigma = sqrt(sigma2);

    tx_bits = randi([0 1], N, 1);
    x = 1 - 2*tx_bits;              % 0 -> +1, 1 -> -1
    y = x + sigma*randn(N,1);

    llr = 2*y/sigma2;
    mag = abs(llr);

    mag_q = round(min(mag, LLR_max));
    mag_q(mag_q > 31) = 31;

    bit_mat = zeros(N, B_mag);
    for b = 1:B_mag
        bit_pos_from_lsb = B_mag - b;
        bit_mat(:,b) = bitget(uint8(mag_q), bit_pos_from_lsb + 1);
    end

    prev_bits = bit_mat(1:end-1,:);
    next_bits = bit_mat(2:end,:);

    P01_all(e,:) = mean(prev_bits == 0 & next_bits == 1, 1);
end

%% Detect transition to white-noise-like region
is_white_noise_like = abs(P01_all - random_ref) < tol;
is_white_noise_like_all = all(is_white_noise_like, 1);

transition_index = find(is_white_noise_like_all, 1, 'first');

if isempty(transition_index)
    transition_x = 2.5;
    useful_bits = 3;
else
    transition_x = transition_index - 0.5;
    useful_bits = transition_index;
end

%% Plot
make_fig(EbN0_list, P01_all, bit_labels, random_ref, transition_x);

%% Summary table
fprintf('\nBit-transition P(0->1) summary:\n');
T = array2table(P01_all, 'VariableNames', ...
    {'MSB4','MSB3','MSB2','LSB1','LSB0'});
T.EbN0_dB = EbN0_list(:);
T = movevars(T, 'EbN0_dB', 'Before', 1);
disp(T);

fprintf('\nTransition analysis:\n');
fprintf('Random-bit reference P(0->1) = %.2f\n', random_ref);
fprintf('Tolerance = %.3f\n', tol);
fprintf('Estimated useful MSB word length = %d bits\n', useful_bits);
fprintf('Useful bits kept: ');
disp(bit_labels(1:useful_bits));

%% ==========================================================
%% Figure function
%% ==========================================================
function make_fig(EbN0_list, P01_all, bit_labels, random_ref, transition_x)

    FONT = 10;
    AXIS_FONT = 8;
    MS = 3.5;
    LW = 1.5;

    colors = [
        0.0000 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.4660 0.6740 0.1880
        0.4940 0.1840 0.5560
    ];

    markers = {'-o','-s','-^','-d'};

    figure(1); clf;
    set(gcf,'Color','w','Position',[100 100 420 300]);

    hold on;

    for e = 1:length(EbN0_list)
        plot(1:size(P01_all,2), P01_all(e,:), markers{e}, ...
            'DisplayName', sprintf('$E_b/N_0=%g$ dB', EbN0_list(e)), ...
            'Color', colors(e,:), ...
            'LineWidth', LW, ...
            'MarkerSize', MS);
    end

    yline(random_ref, 'k--', ...
        'LineWidth', 1.2, ...
        'DisplayName', '$P(0\rightarrow1)=0.25$');

    xline(transition_x, 'r--', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'useful-bit boundary');

    hold off;

    ax = gca;
    ax.FontSize = AXIS_FONT;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    ax.TickLabelInterpreter = 'latex';

    grid on;
    grid minor;
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.1;
    box on;

    xlim([0.7 5.3]);
    ylim([0 0.30]);

    xticks(1:5);
    xticklabels(bit_labels);
    yticks(0:0.05:0.30);

    xlabel('5-bit LLR magnitude bit position','FontSize',FONT);
    ylabel('$P(0 \rightarrow 1)$','FontSize',FONT);

    lgd = legend('show','Location','southeast');
    lgd.FontSize = 7;
    lgd.Box = 'on';

    set(gcf,'Renderer','painters');
    exportgraphics(gcf,'llr_bit_transition.pdf','ContentType','vector');

end