clear codes128 codes256

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% =========================
%  STYLE MAP
%  =========================
sty.Baseline.decoder = 'Baseline';
sty.Baseline.LT      = '-o';
sty.Baseline.color   = [0 0 0];

sty.MSB4.decoder     = 'MSB$_4$';
sty.MSB4.LT          = '-s';
sty.MSB4.color       = [0 0.45 0.74];

sty.MSB3.decoder     = 'MSB$_3$';
sty.MSB3.LT          = '-d';
sty.MSB3.color       = [0.85 0.33 0.10];

sty.PA.decoder       = 'PA';
sty.PA.LT            = '-^';
sty.PA.color         = [0.47 0.67 0.19];

%% =========================
%  LOAD n = 128, k = 116
%  =========================
n = 128;
k = 116;
n_CODES = 0;
codes128 = struct([]);

% ---- PA
DECODER = 'ORBGRAND-MSB3-TIEBREAK2';
code.class = 'CAPOLAR';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = sty.PA.decoder;
code.LT      = sty.PA.LT;
code.color   = sty.PA.color;
codes128(n_CODES).code = code;

% ---- Baseline
DECODER = 'ORBGRAND-BASELINE';
code.class = 'CAPOLAR';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = sty.Baseline.decoder;
code.LT      = sty.Baseline.LT;
code.color   = sty.Baseline.color;
codes128(n_CODES).code = code;

% ---- MSB3
DECODER = 'ORBGRAND-MSB-3';
code.class = 'CAPOLAR';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = sty.MSB3.decoder;
code.LT      = sty.MSB3.LT;
code.color   = sty.MSB3.color;
codes128(n_CODES).code = code;

% ---- MSB4
DECODER = 'ORBGRAND-MSB-4';
code.class = 'CAPOLAR';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = sty.MSB4.decoder;
code.LT      = sty.MSB4.LT;
code.color   = sty.MSB4.color;
codes128(n_CODES).code = code;

%% =========================
%  LOAD n = 256, k = 240
%  =========================
n = 256;
k = 240;
n_CODES = 0;
codes256 = struct([]);

% ---- PA
DECODER = 'ORBGRAND-MSB3-TIEBREAK2';
code.class = 'CRC';
poly = '0xd175';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = sty.PA.decoder;
code.LT      = sty.PA.LT;
code.color   = sty.PA.color;
codes256(n_CODES).code = code;

% ---- Baseline
DECODER = 'ORBGRAND-BASELINE';
code.class = 'CRC';
poly = '0xd175';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = sty.Baseline.decoder;
code.LT      = sty.Baseline.LT;
code.color   = sty.Baseline.color;
codes256(n_CODES).code = code;

% ---- MSB3
DECODER = 'ORBGRAND-MSB-3';
code.class = 'CRC';
poly = '0xd175';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = sty.MSB3.decoder;
code.LT      = sty.MSB3.LT;
code.color   = sty.MSB3.color;
codes256(n_CODES).code = code;

% ---- MSB4
DECODER = 'ORBGRAND-MSB-4';
code.class = 'CRC';
poly = '0xd175';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = sty.MSB4.decoder;
code.LT      = sty.MSB4.LT;
code.color   = sty.MSB4.color;
codes256(n_CODES).code = code;

%% =========================
%  MAKE FIGURE
%  =========================
make_fig_combined(codes128,codes256,1)

function make_fig_combined(codes128,codes256,fig_no)

    FONT = 13;
    AXIS_FONT = 11;
    MS = 6.5;
    LW = 1.4;

    XLIM_ZOOM = [4 7];
    YLIM_FER  = [6e-6 3e-1];
    YTICK_FER = [1e-5 1e-4 1e-3 1e-2 1e-1];

    labels128 = make_labels(codes128);
    labels256 = make_labels(codes256);

    figure(fig_no); clf;
    set(gcf,'Color','w');
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1 1 20.5 12.8]);
    set(gcf,'Renderer','painters');

    t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % ===================== 128 FER =====================
    ax1 = nexttile(1);
    hold(ax1,'on');
    h128 = gobjects(length(codes128),1);
    for ii = 1:length(codes128)
        code = codes128(ii).code;
        h128(ii) = plot(ax1, code.ebn0, code.BLER, code.LT, ...
            'DisplayName', labels128{ii}, ...
            'Color', code.color, ...
            'LineWidth', LW, ...
            'MarkerSize', MS, ...
            'MarkerFaceColor', 'w');
    end
    hold(ax1,'off');
    style_axis_small(ax1,AXIS_FONT,XLIM_ZOOM,YLIM_FER,YTICK_FER);
    xlabel(ax1,'$E_b/N_0$ (dB)','FontSize',FONT);
    ylabel(ax1,'$\mathrm{FER}$','FontSize',FONT);
    title(ax1,'$(128,116)$ CA-Polar','FontSize',FONT);

    % ===================== 128 EG =====================
    ax2 = nexttile(2);
    hold(ax2,'on');
    for ii = 1:length(codes128)
        code = codes128(ii).code;
        plot(ax2, code.ebn0, code.EG, code.LT, ...
            'DisplayName', labels128{ii}, ...
            'Color', code.color, ...
            'LineWidth', LW, ...
            'MarkerSize', MS, ...
            'MarkerFaceColor', 'w');
    end
    hold(ax2,'off');
    style_axis_small(ax2,AXIS_FONT,XLIM_ZOOM,[],[1e0 1e1 1e2 1e3]);
    xlabel(ax2,'$E_b/N_0$ (dB)','FontSize',FONT);
    ylabel(ax2,'Avg. codebook queries','FontSize',FONT);

    % ===================== 256 FER =====================
    ax3 = nexttile(3);
    hold(ax3,'on');
    h256 = gobjects(length(codes256),1);
    for ii = 1:length(codes256)
        code = codes256(ii).code;
        h256(ii) = plot(ax3, code.ebn0, code.BLER, code.LT, ...
            'DisplayName', labels256{ii}, ...
            'Color', code.color, ...
            'LineWidth', LW, ...
            'MarkerSize', MS, ...
            'MarkerFaceColor', 'w');
    end
    hold(ax3,'off');
    style_axis_small(ax3,AXIS_FONT,XLIM_ZOOM,YLIM_FER,YTICK_FER);
    xlabel(ax3,'$E_b/N_0$ (dB)','FontSize',FONT);
    ylabel(ax3,'$\mathrm{FER}$','FontSize',FONT);
    title(ax3,'$(256,240)$ CRC','FontSize',FONT);

    % ===================== 256 EG =====================
    ax4 = nexttile(4);
    hold(ax4,'on');
    for ii = 1:length(codes256)
        code = codes256(ii).code;
        plot(ax4, code.ebn0, code.EG, code.LT, ...
            'DisplayName', labels256{ii}, ...
            'Color', code.color, ...
            'LineWidth', LW, ...
            'MarkerSize', MS, ...
            'MarkerFaceColor', 'w');
    end
    hold(ax4,'off');
    style_axis_small(ax4,AXIS_FONT,XLIM_ZOOM,[],[1e0 1e1 1e2 1e3 1e4 1e5]);
    xlabel(ax4,'$E_b/N_0$ (dB)','FontSize',FONT);
    ylabel(ax4,'Avg. codebook queries','FontSize',FONT);

    % ===================== LEGENDS =====================
    lgd1 = legend(ax2,h128,labels128, ...
        'Interpreter','latex', ...
        'Location','northeast', ...
        'FontSize',8, ...
        'Box','on');
    lgd1.NumColumns = 1;
    lgd1.ItemTokenSize = [14 8];

    lgd2 = legend(ax4,h256,labels256, ...
        'Interpreter','latex', ...
        'Location','northeast', ...
        'FontSize',8, ...
        'Box','on');
    lgd2.NumColumns = 1;
    lgd2.ItemTokenSize = [14 8];

    exportgraphics(gcf, ...
        'orbgrand_sorter_optimization_combined.pdf', ...
        'ContentType','vector', ...
        'BackgroundColor','white');
end

function labels = make_labels(codes)
    labels = cell(1,length(codes));
    for ii = 1:length(codes)
        c = codes(ii).code;
        labels{ii} = c.decoder;
    end
end

function style_axis_small(ax,AXIS_FONT,XLIM,YLIM,YTICKS)
    ax.FontSize = AXIS_FONT;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.Box = 'on';
    ax.YScale = 'log';
    ax.YMinorTick = 'on';

    xlim(ax,XLIM);
    xticks(ax,4:7);
    yticks(ax,YTICKS);

    if ~isempty(YLIM)
        ylim(ax,YLIM);
    end

    grid(ax,'on');
    grid(ax,'minor');
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.10;
end