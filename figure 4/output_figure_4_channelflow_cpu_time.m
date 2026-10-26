%% Plot Memory vs Grid Size for Channel Flow and Wavy Wall (Matrix-free vs Sparse)
clearvars; close all; clc;

folderpath = './'
cd(folderpath)
T = readtable('preconditioner_comparison.xlsx');

%% ---- Style parameters ----
FONT     = 'Arial';
LW       = 3.4;    % line width (also sets marker edge width in MATLAB)
MS       = 24;     % marker size
FS_LABEL = 40;      % axis label font size
FS_TICK  = 44;      % tick label font size
FS_LEG   = 44;      % legend font size

% Your colors
C_MF_CF = [0 0.4470 0.7410];      % blue
%C_LP_CF = [0.8500 0.3250 0.0980]; % orange
C_SP_CF = [0.6350 0.0780 0.1840]; % red
C_LP_CF = [0.1059 0.6863 0.4784];   % teal/green

%% ---- Series masks ----
idx_mf_cf = strcmp(T.MatrixFree_sparseMatrix,'matrix free') & strcmp(T.Channelflow_WavyWall,'Channel flow');
%idx_mf_ww = strcmp(T.MatrixFree_sparseMatrix,'matrix free') & strcmp(T.Channelflow_WavyWall,'Wavy wall');
idx_sp_cf = strcmp(T.MatrixFree_sparseMatrix,'sparse matrix') & strcmp(T.Channelflow_WavyWall,'Channel flow');
idx_lp_cf = strcmp(T.MatrixFree_sparseMatrix,'laplacian') & strcmp(T.Channelflow_WavyWall,'Channel flow');

%% ---- Figure/axes ----
%fig1 = figure('Units','normalized','Position',[0.12 0.12 0.62 0.58],'Color','w');
fig1 = figure('Units','inches','Position',[1 1 16.53 8.7],'Color','w');
ax = axes(fig1); hold(ax,'on'); box(ax,'on');

% Sparse matrix -- channel flow: filled circle, solid line
h_sp_cf = loglog(ax, T.N(idx_sp_cf), T.CPU_sec_(idx_sp_cf), '^-', ...
    'Color', C_SP_CF, 'MarkerFaceColor', C_SP_CF, 'MarkerEdgeColor', C_SP_CF, ...
    'LineWidth', LW, 'MarkerSize', MS);

% Matrix-free -- channel flow: filled circle, solid line
h_mf_cf = loglog(ax, T.N(idx_mf_cf), T.CPU_sec_(idx_mf_cf), 'o-', ...
    'Color', C_MF_CF, 'MarkerFaceColor', C_MF_CF, 'MarkerEdgeColor', C_MF_CF, ...
    'LineWidth', LW, 'MarkerSize', MS);

% Laplacian preconditioner -- channel flow: open square, dashed line
% h_lp_cf = loglog(ax, T.N(idx_lp_cf), T.CPU_sec_(idx_lp_cf), '*--', ...
%     'Color', C_LP_CF, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', C_LP_CF, ...
%     'LineWidth', LW, 'MarkerSize', MS);


%% ---- Axes formatting ----
set(ax, 'XScale', 'log', 'YScale', 'log');

grid(ax, 'on');
ax.GridLineStyle = '-';   ax.GridAlpha = 0.25;
ax.MinorGridLineStyle = ':'; ax.MinorGridAlpha = 0.15;
ax.FontName  = FONT;
ax.FontSize  = FS_TICK;
ax.LineWidth = 1.6;

allN   = T.N;
allCPU = T.CPU_sec_;
xt = 10.^(floor(log10(min(allN))):ceil(log10(max(allN))));
yt = 10.^(floor(log10(min(allCPU))):ceil(log10(max(allCPU))));
set(ax, 'XTick', xt, 'YTick', yt);

xlabel(ax, 'N (grid size)', 'FontSize', FS_LABEL, 'FontWeight', 'bold', 'FontName', FONT);
ylabel(ax, 'CPU time (s)', 'FontSize', FS_LABEL, 'FontWeight', 'bold', 'FontName', FONT);

% legend(ax, [h_sp_cf h_mf_cf h_lp_cf], ...
%     {'Sparse matrix - channel flow', 'Matrix-free - channel flow', ...
%      'Sparse matrix - wavy wall',    'Matrix-free - wavy wall'}, ...
%     'Location', 'southeast', 'FontSize', FS_LEG, 'FontName', FONT, 'Box', 'off');


legend(ax, [h_sp_cf h_mf_cf], ...
    {'Sparse matrix', 'Matrix-free'}, ...
    'Location', 'southeast', 'FontSize', FS_LEG, 'FontName', FONT, 'Box', 'off');

%% ---- Export (raster for slides/preview, vector for LaTeX) ----
exportgraphics(fig1, 'figure 4b.png', 'Resolution', 300);
%exportgraphics(fig1, 'memory_vs_N_paper.pdf', 'ContentType', 'vector');
hold on


%% ---- Series masks ----
idx_mf_cf = strcmp(T.MatrixFree_sparseMatrix,'matrix free') & strcmp(T.Channelflow_WavyWall,'Channel flow');
idx_sp_cf = strcmp(T.MatrixFree_sparseMatrix,'sparse matrix') & strcmp(T.Channelflow_WavyWall,'Channel flow');
%idx_lp_cf = strcmp(T.MatrixFree_sparseMatrix,'laplacian') & strcmp(T.Channelflow_WavyWall,'Channel flow');
 
close all;