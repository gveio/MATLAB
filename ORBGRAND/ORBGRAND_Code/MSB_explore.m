clear codes

% --- LaTeX style ---
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

n = 256;
k = 240;
n_CODES = 0;

code.class = 'CRC';
poly='0xd175';
% -------------------------
% Load baseline
% -------------------------
DECODER = 'ORBGRAND-BASELINE';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = 'Baseline';
code.LT = '-d';
code.color = 'k';
codes(n_CODES).code = code;

% -------------------------
% Load MSB2
% -------------------------
DECODER = 'ORBGRAND-MSB-2';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = 'MSB$_2$';
code.LT = '-o';
code.color = [0.8500 0.3250 0.0980];
codes(n_CODES).code = code;

% -------------------------
% Load MSB3
% -------------------------
DECODER = 'ORBGRAND-MSB-3';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = 'MSB$_3$';
code.LT = '-s';
code.color = [0 0.4470 0.7410];
codes(n_CODES).code = code;

% -------------------------
% Load MSB4
% -------------------------
DECODER = 'ORBGRAND-MSB-4';
n_CODES = n_CODES + 1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = 'MSB$_4$';
code.LT = '-^';
code.color = [0.4660 0.6740 0.1880];
codes(n_CODES).code = code;

% -------------------------
% Make figure
% -------------------------
make_fig(codes,1)

function make_fig(codes,fig_no)

    FONT = 10;
    AXIS_FONT = 8;
    MS = 3.5;
    LW = 1.5;

    XLIM_ZOOM = [4 7];
    YLIM_ZOOM = [1e-5 1];

    num_codes = length(codes);

    figure(fig_no); clf;
    set(gcf,'Color','w','Position',[100 100 420 300]); % smaller figure

    hold on
    for ii = 1:num_codes
        plot(codes(ii).code.ebn0, ...
             codes(ii).code.BLER, ...
             codes(ii).code.LT, ...
             'DisplayName', codes(ii).code.decoder, ...
             'Color', codes(ii).code.color, ...
             'LineWidth', LW, ...
             'MarkerSize', MS);
    end
    hold off

    ax = gca;
    ax.FontSize = AXIS_FONT;
    ax.LineWidth = 1;
    ax.TickDir = 'out';

    set(gca,'YScale','log');
    xlim(XLIM_ZOOM);
    ylim(YLIM_ZOOM);
    xticks(4:7);

    xlabel('$E_b/N_0$ (dB)','FontSize',FONT);
    ylabel('$\mathrm{FER}$','FontSize',FONT);

  lgd = legend('show','Location','northeast');
lgd.FontSize = 8;
lgd.Box = 'on';

   grid on
grid minor
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.1;
    box on

    % ---------- Export ----------
    set(gcf,'Renderer','painters');
    exportgraphics(gcf,'msb_precision_fer.pdf','ContentType','vector');
    % ----------------------------

end