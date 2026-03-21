 % ==========================================================
% plot_tie_probability_figures.m
%
% Creates paper-ready IEEE SiPS single-column PDF figures
% from precomputed tie-probability .mat files.
%
% Required input files:
%   ../ORBGRAND_Code/TIE_STATS_CAPOLAR_128_116.mat
%   ../ORBGRAND_Code/TIE_STATS_CRC_0xd175_256_240.mat
%
% Output figures:
%   ../ORBGRAND_Code/tie_probability_stage_subplots.pdf
%   ../ORBGRAND_Code/tie_probability_normalized_combined.pdf
% ==========================================================

clear; clc; close all;

%% ----------------------------------------------------------
% Global plotting defaults
% -----------------------------------------------------------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ----------------------------------------------------------
% Paths
% -----------------------------------------------------------
results_dir = '../ORBGRAND_Code';

file128 = fullfile(results_dir,'TIE_STATS_CAPOLAR_128_116.mat');
file256 = fullfile(results_dir,'TIE_STATS_CRC_0xd175_256_240.mat');

if ~exist(file128,'file')
    error('Missing file: %s', file128);
end
if ~exist(file256,'file')
    error('Missing file: %s', file256);
end

S128 = load(file128);
S256 = load(file256);

result128 = S128.result;
result256 = S256.result;

%% ----------------------------------------------------------
% Basic consistency checks
% -----------------------------------------------------------
if ~isequal(result128.ebn0, result256.ebn0)
    error('The Eb/N0 vectors of the two result files do not match.');
end

ebn0   = result128.ebn0(:).';
numSNR = numel(ebn0);

%% ----------------------------------------------------------
% Plot styling
% -----------------------------------------------------------
colors   = lines(numSNR);

lw_main  = 0.9;
ms_main  = 4.0;
ax_fs    = 8;
lab_fs   = 9;
leg_fs   = 7;
title_fs = 9;

% IEEE single-column size
figW  = 3.5;   % inches
figH1 = 2.25;  % figure 1
figH2 = 2.35;  % figure 2

%% ----------------------------------------------------------
% FIGURE 1:
% Two subplots with raw stage axis
% (a) 128, (b) 256
% -----------------------------------------------------------
f1 = figure('Units','inches', ...
            'Position',[1 1 figW figH1], ...
            'Color','w', ...
            'PaperPositionMode','auto');

tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ----- (a) 128
ax1 = nexttile;
hold(ax1,'on');
for ii = 1:numSNR
    plot(ax1, 1:result128.n_stages, result128.tie_prob(:,ii), '-o', ...
        'Color', colors(ii,:), ...
        'LineWidth', lw_main, ...
        'MarkerSize', ms_main, ...
        'DisplayName', sprintf('$E_b/N_0=%d$ dB', ebn0(ii)));
end
hold(ax1,'off');
grid(ax1,'on');
box(ax1,'on');
ax1.FontSize = ax_fs;
ax1.LineWidth = 0.8;
xlim(ax1,[1 result128.n_stages]);
xticks(ax1,1:result128.n_stages);
ylim(ax1,[0 1]);
xlabel(ax1,'Bitonic stage','FontSize',lab_fs);
ylabel(ax1,'Tie probability','FontSize',lab_fs);
title(ax1,'(a) n = 128','FontSize',title_fs);

% ----- (b) 256
ax2 = nexttile;
hold(ax2,'on');
for ii = 1:numSNR
    plot(ax2, 1:result256.n_stages, result256.tie_prob(:,ii), '-o', ...
        'Color', colors(ii,:), ...
        'LineWidth', lw_main, ...
        'MarkerSize', ms_main, ...
        'DisplayName', sprintf('$E_b/N_0=%d$ dB', ebn0(ii)));
end
hold(ax2,'off');
grid(ax2,'on');
box(ax2,'on');
ax2.FontSize = ax_fs;
ax2.LineWidth = 0.8;
xlim(ax2,[1 result256.n_stages]);
xticks(ax2,1:result256.n_stages);
ylim(ax2,[0 1]);
xlabel(ax2,'Bitonic stage','FontSize',lab_fs);
ylabel(ax2,'Tie probability','FontSize',lab_fs);
title(ax2,'(b) n = 256','FontSize',title_fs);

lg1 = legend(ax2,'show', ...
    'Location','southoutside', ...
    'NumColumns',2, ...
    'FontSize',leg_fs, ...
    'Box','off');
lg1.Layout.Tile = 'south';

exportgraphics(f1, fullfile(results_dir,'tie_probability_stage_subplots.pdf'), ...
    'ContentType','vector');

fprintf('Saved: %s\n', fullfile(results_dir,'tie_probability_stage_subplots.pdf'));

%% ----------------------------------------------------------
% FIGURE 2:
% Combined plot with normalized stage
% color  = Eb/N0
% solid  = 128
% dashed = 256
% -----------------------------------------------------------
f2 = figure('Units','inches', ...
            'Position',[1 1 figW figH2], ...
            'Color','w', ...
            'PaperPositionMode','auto');

% Main axes
ax = axes(f2);
hold(ax,'on');

xnorm128 = (1:result128.n_stages) / result128.n_stages;
xnorm256 = (1:result256.n_stages) / result256.n_stages;

for ii = 1:numSNR
    % n = 128 -> solid with circle
    plot(ax, xnorm128, result128.tie_prob(:,ii), '-o', ...
        'Color', colors(ii,:), ...
        'LineWidth', lw_main, ...
        'MarkerSize', ms_main);

    % n = 256 -> dashed with square
    plot(ax, xnorm256, result256.tie_prob(:,ii), '--s', ...
        'Color', colors(ii,:), ...
        'LineWidth', lw_main, ...
        'MarkerSize', ms_main-0.2);
end

grid(ax,'on');
box(ax,'on');
ax.FontSize = ax_fs;
ax.LineWidth = 0.8;

xlim(ax,[min([xnorm128(:); xnorm256(:)]) 1.0]);
xticks(ax,[0.25 0.50 0.75 1.00]);
ylim(ax,[0 1]);

xlabel(ax,'Normalized stage','FontSize',lab_fs);
ylabel(ax,'Tie probability','FontSize',lab_fs);

% ----------------------------------------------------------
% Bottom legend: colors only (no markers)
% ----------------------------------------------------------
hColor = gobjects(numSNR,1);
for ii = 1:numSNR
    hColor(ii) = plot(ax, nan, nan, '-', ...
        'Color', colors(ii,:), ...
        'LineWidth', lw_main);
end

lg1 = legend(ax, hColor, ...
    compose('$E_b/N_0=%d$ dB', ebn0), ...
    'Location','southoutside', ...
    'NumColumns',2, ...
    'FontSize',leg_fs, ...
    'Box','off');

% ----------------------------------------------------------
% Second invisible axes only for line-style + marker legend
% ----------------------------------------------------------
ax2 = axes('Position', ax.Position, ...
           'Color', 'none', ...
           'XColor', 'none', ...
           'YColor', 'none', ...
           'HitTest', 'off', ...
           'HandleVisibility', 'off');
hold(ax2,'on');

h1 = plot(ax2, nan, nan, '-o', ...
    'Color','k', ...
    'LineWidth',lw_main, ...
    'MarkerSize',ms_main);

h2 = plot(ax2, nan, nan, '--s', ...
    'Color','k', ...
    'LineWidth',lw_main, ...
    'MarkerSize',ms_main-0.2);

lg2 = legend(ax2, [h1 h2], {'$n=128$','$n=256$'}, ...
    'Location','northwest', ...
    'FontSize',leg_fs, ...
    'Box','on');

hold(ax2,'off');

exportgraphics(f2, fullfile(results_dir,'tie_probability_normalized_combined.pdf'), ...
    'ContentType','vector');

fprintf('Saved: %s\n', fullfile(results_dir,'tie_probability_normalized_combined.pdf'));