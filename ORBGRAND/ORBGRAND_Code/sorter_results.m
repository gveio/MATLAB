%% ORBGRAND – Final Hardware Results (Thesis Version)
clear; clc; close all;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesFontSize',12);

Architecture = {'Baseline ORBGRAND'; 'PA-ORBGRAND'};

N = [256; 256];
K = [240; 240];

%% Sorter Resources

Sorter_LUTs = [93440; 76139];
Sorter_FFs  = [21758; 18526];

%% Power Results

% Full-chip power (PS+PL)
DynPower_total   = [3.567; 2.798];
TotalPower_total = [4.265; 3.490];

% PL-only (decoder wrapper)
PL_Dynamic = [0.929; 0.158];
PL_Total   = [1.527; 0.752];

% Sorter-only dynamic
DynPower_Sorter = [0.869; 0.087];

%% Percentage Reductions

LUT_red  = 100*(Sorter_LUTs(1)-Sorter_LUTs(2))/Sorter_LUTs(1);
FF_red   = 100*(Sorter_FFs(1)-Sorter_FFs(2))/Sorter_FFs(1);
PL_red   = 100*(PL_Total(1)-PL_Total(2))/PL_Total(1);
SorterP_red = 100*(DynPower_Sorter(1)-DynPower_Sorter(2))/DynPower_Sorter(1);

fprintf('\nHardware Savings:\n');
fprintf('Sorter LUT reduction: %.2f %%\n', LUT_red);
fprintf('Sorter FF reduction: %.2f %%\n', FF_red);
fprintf('PL Total Power reduction: %.2f %%\n', PL_red);
fprintf('Sorter Dynamic reduction: %.2f %%\n', SorterP_red);

%% Sorter Area Summary Plot

figure;
bar([Sorter_LUTs Sorter_FFs],'grouped')
grid on
ylabel('Resources')
legend('LUTs','FFs','Location','northoutside','Orientation','horizontal')
title('Bitonic Sorter Resource Utilization')
xticklabels(Architecture)
xtickangle(20)

text(2,Sorter_LUTs(2),sprintf('-%.1f%% LUT',LUT_red),'VerticalAlignment','bottom')
text(2,Sorter_FFs(2),sprintf('-%.1f%% FF',FF_red),'VerticalAlignment','bottom')

exportgraphics(gcf,'sorter_area_summary.png','ContentType','vector');

%% Power Summary Plot

figure;
bar([PL_Total DynPower_Sorter],'grouped')
grid on
ylabel('Power (W)')
legend('PL Total','Sorter Dynamic','Location','northoutside','Orientation','horizontal')
title('Power Reduction Summary (115 MHz)')
xticklabels(Architecture)
xtickangle(20)

text(2,PL_Total(2),sprintf('-%.1f%% PL',PL_red),'VerticalAlignment','bottom')
text(2,DynPower_Sorter(2),sprintf('-%.1f%% Sorter',SorterP_red),'VerticalAlignment','bottom')

exportgraphics(gcf,'power_summary.png','ContentType','vector');

%% Normalized Hardware Comparison

Norm_Area  = Sorter_LUTs ./ Sorter_LUTs(1);
Norm_Power = PL_Total ./ PL_Total(1);

figure;
bar([Norm_Area Norm_Power],'grouped')
grid on
ylabel('Normalized (Baseline = 1)')
legend('Sorter LUTs','PL Total Power','Location','northoutside','Orientation','horizontal')
title('Normalized Hardware Reduction')
xticklabels(Architecture)
xtickangle(20)

exportgraphics(gcf,'normalized_summary.png','ContentType','vector');