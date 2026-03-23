% ORBGRAND area pie (wide layout for manual annotations)

clear; clc; close all;

% ---------------- Data ----------------
area_sorter        = 24756.694;
area_error_gen     = 1189.189;
area_pattern_gen   = 345.531;
area_total         = 28026.726;

area_error_checker = area_total - (area_sorter + area_error_gen + area_pattern_gen);

vals = [area_error_checker, area_sorter, area_pattern_gen, area_error_gen];

% ---------------- Figure ----------------
fig = figure('Color','white');

% Make figure WIDE (important)
set(fig, 'Units', 'centimeters');
set(fig, 'Position', [2 2 16 10]);  % <-- wide layout

ax = axes(fig);

% Shrink pie inside the figure (leave margins for labels)
ax.Position = [0.25 0.15 0.5 0.7];  
% [left bottom width height]

hold(ax,'on');

% Draw pie
p = pie(ax, vals);

% ---------------- Colors ----------------
% Ensure correct order of patches
patchHandles = flipud(findobj(p, 'Type', 'Patch'));

colors = [
    0.88 0.88 0.88;   % Error Checker
    0.85 0.75 0.45;   % Pattern Generator (soft yellow)
    0.35 0.55 0.75;   % Bitonic Sorter (muted blue)
    0.75 0.45 0.45;   % Error Generator (soft red)
];

for k = 1:numel(patchHandles)
    patchHandles(k).FaceColor = colors(k,:);
    patchHandles(k).EdgeColor = [0.4 0.4 0.4];
    patchHandles(k).LineWidth = 0.6;
end

% ---------------- Remove labels ----------------
textHandles = findobj(p, 'Type', 'Text');
delete(textHandles);

% ---------------- Clean layout ----------------
axis(ax,'equal');
axis(ax,'off');

% ---------------- Export ----------------
exportgraphics(fig, 'orbgrand_area.png', 'ContentType', 'vector');