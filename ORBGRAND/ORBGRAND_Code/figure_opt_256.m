clear codes

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');


n=256;
k=240;
n_CODES = 0;

% ORBGRAND
DECODER = 'ORBGRAND-MSB-3';

code.class = 'CRC';
poly='0xd175';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = [DECODER ' ' poly];
code.LT = '-^';
code.color = 'r';
codes(n_CODES).code = code; 


DECODER = 'ORBGRAND-MSB3-TIE-FLAG';

code.class = 'CRC';
poly='0xd175';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = [DECODER ' ' poly];
code.LT = '-v';
code.color = 'r';
codes(n_CODES).code = code; 

DECODER = 'ORBGRAND-MSB3-TIEBREAK1';

code.class = 'CRC';
poly='0xd175';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = [DECODER ' ' poly];
code.LT = '-O';
code.color = 'r';
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
    lgd.Box = 'on';
lgd.Interpreter = 'tex';   % ← NOT latex

drawnow;  % force layout calculation

lgd.Units = 'normalized';
pos = lgd.Position;

pos(1) = 0.5 - pos(3)/2;   % EXACT horizontal center
pos(2) = 0.93;             % vertical height (adjust slightly if needed)

lgd.Position = pos;


    % ===================== EXPORT =====================
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1 1 22 7]);        % wide figure
    set(gcf,'Renderer','painters');

    drawnow;                                  % sync graphics

    set(gcf,'Units','normalized');
    set(gcf,'Position',[0 0 1 1]);            % FULL SCREEN

exportgraphics(gcf, ...
    'orbgrand_sorter_optimization_256.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');

end
