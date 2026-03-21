% Sample plotted output from simulations with [128,116] codes and ORBGRAND (soft detection)
% Silulation of C implementation
clear codes

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

n=128;
k=116;
n_CODES = 0;

% ORBGRAND (soft detection)

DECODER = 'ORBGRAND_C';

% code.class = 'BCH';
% n_CODES=n_CODES+1;
% filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n-1) '_' num2str(k-3) '_1.mat'];
% load(filename,'code');
% code.decoder = DECODER;
% code.LT = '-d';
% code.color = 'g';
% codes(n_CODES).code = code; 

code.class = 'CAPOLAR';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color = 'b';
codes(n_CODES).code = code; 

code.class = 'PAC';
conv_code = [1 1 1 0 1 1 0 1];
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(bin2dec(num2str(conv_code))) '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-s';
code.color = [0.1250, 0.6940, 0.9290];
codes(n_CODES).code = code; 

code.class = 'CRC';
poly='0x8f3';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = [DECODER ' ' poly];
code.LT = '-o';
code.color = 'r';
codes(n_CODES).code = code; 

code.class = 'RLC';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-x';
code.color = 'm';
codes(n_CODES).code = code; 

DECODER = 'ORBGRAND_M';

% code.class = 'BCH';
% n_CODES=n_CODES+1;
% filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n-1) '_' num2str(k-3) '_1.mat'];
% load(filename,'code');
% code.decoder = DECODER;
% code.LT = '-.d';
% code.color = 'g';
% codes(n_CODES).code = code; 

code.class = 'CAPOLAR';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-.d';
code.color = 'b';
codes(n_CODES).code = code; 

code.class = 'PAC';
conv_code = [1 1 1 0 1 1 0 1];
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(bin2dec(num2str(conv_code))) '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-.s';
code.color = [0.1250, 0.6940, 0.9290];
codes(n_CODES).code = code; 

code.class = 'CRC';
poly='0x8f3';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = [DECODER ' ' poly];
code.LT = '-.o';
code.color = 'r';
codes(n_CODES).code = code; 

code.class = 'RLC';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-.x';
code.color = 'm';
codes(n_CODES).code = code; 

make_fig(codes,1)

function make_fig(codes,fig_no)

    FONT = 14;
    MS   = 7;
    LW   = 1.8;
    num_codes = length(codes);
    labels = cell(1,num_codes);
    for ii = 1:num_codes
         c = codes(ii).code;
         labels{ii} = sprintf('%s, %s [%d,%d], R=%.2f', ...
         c.decoder, c.class, c.n, c.k, c.k/c.n);
    end

    figure(fig_no); clf;

    % ---- Layout ----
    t = tiledlayout(1,3, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    % ===================== BLER =====================
    nexttile
    hold on
    for ii = 1:num_codes
        code = codes(ii).code;
        plot(code.ebn0, code.BLER, code.LT, ...
        'DisplayName', labels{ii}, ...
        'Color', code.color, ...
        'LineWidth', LW, 'MarkerSize', MS);

    end
    hold off
    grid on
    ylabel('BLER')
    xlabel('$E_b/N_0$ (dB)')
    set(gca,'YScale','log','FontSize',FONT)

    % ===================== BER =====================
    nexttile
    hold on
    for ii = 1:num_codes
        code = codes(ii).code;
         plot(code.ebn0, code.BER, code.LT, ...
        'DisplayName', labels{ii}, ...
        'Color', code.color, ...
        'LineWidth', LW, 'MarkerSize', MS);

    end
    hold off
    grid on
    ylabel('BER')
    xlabel('$E_b/N_0$ (dB)')
    set(gca,'YScale','log','FontSize',FONT)

    % ===================== EG =====================
    nexttile
    hold on
    for ii = 1:num_codes
        code = codes(ii).code;
        plot(code.ebn0, code.EG, code.LT, ...
        'DisplayName', labels{ii}, ...
        'Color', code.color, ...
        'LineWidth', LW, 'MarkerSize', MS);

    end
    hold off
    grid on
    ylabel('Avg. code-book queries')
    xlabel('$E_b/N_0$ (dB)')
    set(gca,'YScale','log','FontSize',FONT)

    % ===================== SHARED LEGEND =====================
    lgd = legend('show');
    lgd.Layout.Tile = 'north';
    lgd.FontSize = 10;

    % ===================== EXPORT =====================
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1 1 22 7]);        % wide figure
    set(gcf,'Renderer','painters');

    drawnow;                                  % sync graphics

    set(gcf,'Units','normalized');
    set(gcf,'Position',[0 0 1 1]);            % FULL SCREEN

exportgraphics(gcf, ...
    'orbgrand_cmex_vs_matlab.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');

end
