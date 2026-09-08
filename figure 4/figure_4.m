%% Plot Memory vs Grid Size for Channel Flow and Wavy Wall (Matrix-free vs Sparse)
clearvars; close all; clc;

folderpath = './'
cd(folderpath)
T = readtable('job id 1.xlsx');

%% ---- Style parameters ----
FONT     = 'Arial';
LW       = 2.4;    % line width (also sets marker edge width in MATLAB)
MS       = 24;     % marker size
FS_LABEL = 40;      % axis label font size
FS_TICK  = 44;      % tick label font size
FS_LEG   = 30;      % legend font size

% Your colors
C_MF_CF = [0 0.4470 0.7410];      % blue
C_MF_WW = [0.8500 0.3250 0.0980]; % orange
C_SP_CF = [0.6350 0.0780 0.1840]; % red
C_SP_WW = [0.1059 0.6863 0.4784];   % teal/green


%% ---- Series masks ----
idx_mf_cf = strcmp(T.MatrixFree_sparseMatrix,'matrix free') & strcmp(T.Channelflow_WavyWall,'Channel flow');
idx_mf_ww = strcmp(T.MatrixFree_sparseMatrix,'matrix free') & strcmp(T.Channelflow_WavyWall,'Wavy wall');
idx_sp_cf = strcmp(T.MatrixFree_sparseMatrix,'sparse matrix') & strcmp(T.Channelflow_WavyWall,'Channel flow');
idx_sp_ww = strcmp(T.MatrixFree_sparseMatrix,'sparse matrix') & strcmp(T.Channelflow_WavyWall,'Wavy wall');

%% ---- Figure/axes ----
fig1 = figure('Units','normalized','Position',[0.12 0.12 0.62 0.58],'Color','w');
ax = axes(fig1); hold(ax,'on'); box(ax,'on');

% Sparse matrix -- channel flow: filled circle, solid line
h_sp_cf = loglog(ax, T.N(idx_sp_cf), T.memory_GB_(idx_sp_cf), '^-', ...
    'Color', C_SP_CF, 'MarkerFaceColor', C_SP_CF, 'MarkerEdgeColor', C_SP_CF, ...
    'LineWidth', LW, 'MarkerSize', MS);

% Matrix-free -- channel flow: filled circle, solid line
h_mf_cf = loglog(ax, T.N(idx_mf_cf), T.memory_GB_(idx_mf_cf), 'o-', ...
    'Color', C_MF_CF, 'MarkerFaceColor', C_MF_CF, 'MarkerEdgeColor', C_MF_CF, ...
    'LineWidth', LW, 'MarkerSize', MS);

% Sparse matrix -- wavy wall: open square, dashed line
h_sp_ww = loglog(ax, T.N(idx_sp_ww), T.memory_GB_(idx_sp_ww), '*--', ...
    'Color', C_SP_WW, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', C_SP_WW, ...
    'LineWidth', LW, 'MarkerSize', MS);

% Matrix-free -- wavy wall: open square, dashed line
h_mf_ww = loglog(ax, T.N(idx_mf_ww), T.memory_GB_(idx_mf_ww), 's--', ...
    'Color', C_MF_WW, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', C_MF_WW, ...
    'LineWidth', LW, 'MarkerSize', MS);

%% ---- Axes formatting ----
set(ax, 'XScale', 'log', 'YScale', 'log');
grid(ax, 'on');
ax.GridLineStyle = '-';   ax.GridAlpha = 0.25;
ax.MinorGridLineStyle = ':'; ax.MinorGridAlpha = 0.15;
ax.FontName  = FONT;
ax.FontSize  = FS_TICK;
ax.LineWidth = 1.6;

allN   = T.N;
allMem = T.memory_GB_;
xt = 10.^(floor(log10(min(allN))):ceil(log10(max(allN))));
yt = 10.^(floor(log10(min(allMem))):ceil(log10(max(allMem))));
set(ax, 'XTick', xt, 'YTick', yt);

xlabel(ax, 'N (grid size)', 'FontSize', FS_LABEL, 'FontWeight', 'bold', 'FontName', FONT);
ylabel(ax, 'Memory (GB)', 'FontSize', FS_LABEL, 'FontWeight', 'bold', 'FontName', FONT);

legend(ax, [h_sp_cf h_mf_cf h_sp_ww h_mf_ww], ...
    {'Sparse matrix - channel flow', 'Matrix-free - channel flow', ...
     'Sparse matrix - wavy wall',    'Matrix-free - wavy wall'}, ...
    'Location', 'southeast', 'FontSize', FS_LEG, 'FontName', FONT, 'Box', 'off');

%% ---- Export (raster for slides/preview, vector for LaTeX) ----
exportgraphics(fig1, 'figure 4.png', 'Resolution', 300);
%exportgraphics(fig1, 'memory_vs_N_paper.pdf', 'ContentType', 'vector');
hold on


%% ---- Series masks ----
idx_mf_cf = strcmp(T.MatrixFree_sparseMatrix,'matrix free') & strcmp(T.Channelflow_WavyWall,'Channel flow');
idx_mf_ww = strcmp(T.MatrixFree_sparseMatrix,'matrix free') & strcmp(T.Channelflow_WavyWall,'Wavy wall');
idx_sp_cf = strcmp(T.MatrixFree_sparseMatrix,'sparse matrix') & strcmp(T.Channelflow_WavyWall,'Channel flow');
idx_sp_ww = strcmp(T.MatrixFree_sparseMatrix,'sparse matrix') & strcmp(T.Channelflow_WavyWall,'Wavy wall');
 
close all;