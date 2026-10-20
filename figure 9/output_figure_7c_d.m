%% ========================================================================
%  Input-Output Postprocessing - Symmetric Wavy Wall
%  Sparse matrix
%
%  Response: U, V, W
%  Forcing : X, Y, Z
%% ========================================================================
clear all;
clearvars;
close all;
clc;
%% ========================= USER INPUTS ==================================
folderpath = './'
results_file = ...
    'stability_results_Re500_kz10_156x116.mat';
cd(folderpath)
%% ========================= PARAMETERS ===================================
Ny = 156;
Nx = 116;
epsilon = 0.12;
Re = 500;

y0 = 1;

A1 = epsilon;
A2 = epsilon;

Lx = 2.4;
%% ========================= PLOT CONTROLS ================================

show_velocity = true;

% Font sizes -- same as resolvent postprocessing
axis_tick_fontsize     = 60;
axis_label_fontsize    = 64;
colorbar_tick_fontsize = 60;

% Solid-region display
solid_color = [0.30 0.30 0.30];

solid_boundary_color = 'k';
solid_boundary_linewidth = 2.0;

% Same figure size as resolvent postprocessing
figure_position = [100 100 1350 1000];

% Number of contour levels
n_contours = 40;

% Output resolution
output_resolution = 300;

%% ========================= LOAD RESULTS =================================
S = load(fullfile(folderpath,results_file));

x = S.x;
y = S.y;

% Response modes
U_hat = S.U_hat;
V_hat = S.V_hat;
W_hat = S.W_hat;

% Forcing modes
U2_hat = S.U2_hat;
V2_hat = S.V2_hat;
W2_hat = S.W2_hat;

% Other result quantities
kz_list = S.kz_list;
c_list  = S.c_list;

if isfield(S,'result_sigma')
    result_sigma = S.result_sigma;
end
%% ========================= KZ INDEX =====================================
tokens = regexp(results_file,'kz(\d+)','tokens');

if isempty(tokens)
    error('Could not determine kz index from filename.');
end

kz_idx = str2double(tokens{1});

kz = kz_list(kz_idx);

fprintf('kz index = %d\n',kz_idx);
fprintf('kz       = %.6g\n',kz);
%% ========================= LOAD SOLID MASK ==============================
filename_mask = sprintf( ...
    'mask_smooth_wavywall_Re%d_%dx%d_epsilon_%.2f.mat', ...
    Re,Ny,Nx,epsilon);

M = load(fullfile(folderpath,filename_mask));
mask_smooth = M.mask_smooth;
solid_mask = (mask_smooth == Re);
%% ========================= CHECK DIMENSIONS =============================
if size(U_hat,1) ~= Ny || size(U_hat,2) ~= Nx
    error('U_hat size does not match Ny x Nx.');
end

if ~isequal(size(solid_mask),[Ny Nx])
    error('solid_mask must have size Ny x Nx.');
end


%% ========================= OUTPUT FOLDER ================================
snapdir = 'output_figure_6c_d';
if ~exist(snapdir,'dir')
    mkdir(snapdir);
end
%% ========================= COORDINATES ==================================
x = x(:).';
y = y(:);
%% ========================= WAVY-WALL CURVES =============================
% Same x-grid as contour plots
x_wall = x(:).';

y1_wall = -y0 - A1*sin(2*pi/Lx*x_wall);
y2_wall =  y0 + A2*sin(2*pi/Lx*x_wall);
%% ========================= DATA ARRAYS ==================================
% Response
resp_cells = {U_hat,V_hat,W_hat};
resp_names = {'U','V','W'};

% Forcing
forc_cells = {U2_hat,V2_hat,W2_hat};
forc_names = {'X','Y','Z'};
%% ========================= C INDICES ====================================
n_modes = size(U_hat,3);

idx_list = 1:n_modes;

%% ========================= COMMON COLOR LIMITS ==========================
resp_clim = zeros(1,3);
forc_clim = zeros(1,3);

for k = 1:3

    %% Response
    for i = idx_list

        tmp = resp_cells{k}(:,:,i);

        tmp(solid_mask) = NaN;

        resp_clim(k) = max( ...
            resp_clim(k), ...
            max(abs(tmp(:)),[],'omitnan'));
    end

    %% Forcing

    for i = idx_list

        tmp = forc_cells{k}(:,:,i);

        tmp(solid_mask) = NaN;

        forc_clim(k) = max( ...
            forc_clim(k), ...
            max(abs(tmp(:)),[],'omitnan'));
    end

end

fprintf('\nCommon response color limits:\n');
disp(resp_clim)

fprintf('Common forcing color limits:\n');
disp(forc_clim)
%% ========================================================================
%  RESPONSE PLOTS
%% ========================================================================
if show_velocity

    for i = 1 %idx_list

        fprintf('Plotting c index %d / %d\n',i,n_modes);


        for k = 1:1

            %% Get field
            field_plot = resp_cells{k}(:,:,i);

            % Hide solid region
            field_plot(solid_mask) = NaN;
            %% ------------------------------------------------------------
            % FIGURE
            % EXACT SAME STYLE AS postprocess_resolvent_wavywall
            % ------------------------------------------------------------

            fig = figure( ...
                'Visible','off', ...
                'Units','pixels', ...
                'Position',figure_position);
            %% Contour
            contourf( ...
                x, ...
                y, ...
                field_plot, ...
                n_contours);

            axis xy;
            hold on;

            clim_val = resp_clim(k);

            ax = gca;

            clim(ax,[-clim_val,clim_val]);


            %% Colorbar
            %
            % clim is intentionally set BEFORE colorbar.
            cb = colorbar;

            cb.Ticks = linspace( ...
                -clim_val, ...
                 clim_val, ...
                 7);

            cb.TickLabels = strip_zeros( ...
                compose('%.3f',cb.Ticks));


            %% Colormap
            colormap(ax,bluewhitered(256));
            %% Domain limits BEFORE wall fill
            xlim(ax,[0 Lx]);
            ylim(ax,[min(y) max(y)]);

            %% Solid wavy-wall patches

            fill_patch_wavywall( ...
                x_wall, ...
                y1_wall, ...
                y2_wall, ...
                solid_color);


            %% Wall boundaries
            plot( ...
                x_wall, ...
                y1_wall, ...
                solid_boundary_color, ...
                'LineWidth',solid_boundary_linewidth);

            plot( ...
                x_wall, ...
                y2_wall, ...
                solid_boundary_color, ...
                'LineWidth',solid_boundary_linewidth);


            %% Domain limits
            %
            % Same as resolvent postprocessing.

            xlim(ax,[0 Lx]);
            ylim(ax,[min(y) max(y)]);


            %% Ticks
            %
            % EXACT SAME layout as resolvent plots.

            ax.XTick = linspace( ...
                0, ...
                Lx, ...
                4);

            ax.YTick = linspace( ...
                min(y), ...
                max(y), ...
                5);


            %% Tick labels -- one decimal

            ax.XTickLabel = strip_zeros( ...
                compose('%.1f',ax.XTick));

            ax.YTickLabel = strip_zeros( ...
                compose('%.1f',ax.YTick));


            %% Fonts

            set( ...
                ax, ...
                'FontSize',axis_tick_fontsize);

            set( ...
                cb, ...
                'FontSize',colorbar_tick_fontsize);


            xlabel( ...
                ax, ...
                'x', ...
                'FontSize',axis_label_fontsize);

            ylabel( ...
                ax, ...
                'y', ...
                'FontSize',axis_label_fontsize);


            %% ------------------------------------------------------------
            % AXES CENTERING
            %
            % DO NOT CHANGE THIS.
            %
            % This is copied directly from the working resolvent
            % postprocessor.
            % ------------------------------------------------------------

%             set( ...
%                 ax, ...
%                 'Position',[0.22 0.20 0.46 0.54]);

            set(ax,'Position',[0.22 0.28 0.46 0.50]);
            %% Physical aspect ratio

            daspect([1 1 1]);

            set(fig,'PaperPositionMode','auto');


            %% ========================= SAVE CSV ==========================

            csv_filename = sprintf( ...
                '%s_response_c%02d_kz%g.csv', ...
                resp_names{k}, ...
                i, ...
                kz);

            writematrix( ...
                field_plot, ...
                fullfile(snapdir,csv_filename));


            %% ========================= SAVE PNG ==========================

            png_filename = sprintf( ...
                '%s_response_c%02d_kz%g.png', ...
                resp_names{k}, ...
                i, ...
                kz);

            exportgraphics( ...
                fig, ...
                fullfile(snapdir,png_filename), ...
                'Resolution',output_resolution, ...
                'BackgroundColor','white');

            close(fig);

        end

        %% =================================================================
        %  FORCING PLOTS
        %% =================================================================

        for k = 1:1

            %% Get field

            field_plot = forc_cells{k}(:,:,i);

            field_plot(solid_mask) = NaN;


            %% Figure
            %
            % SAME as response / resolvent plot.

            fig = figure( ...
                'Visible','off', ...
                'Units','pixels', ...
                'Position',figure_position);


            %% Contour

            contourf( ...
                x, ...
                y, ...
                field_plot, ...
                n_contours);

            axis xy;
            hold on;


            %% Color limits

            clim_val = forc_clim(k);

            ax = gca;

            clim(ax,[-clim_val,clim_val]);


            %% Colorbar
            cb = colorbar;

            cb.Ticks = linspace( ...
                -clim_val, ...
                 clim_val, ...
                 7);

            cb.TickLabels = strip_zeros( ...
                compose('%.3f',cb.Ticks));


            %% Colormap
            colormap(ax,bluewhitered(256));


            %% Domain limits BEFORE wall fill

            xlim(ax,[0 Lx]);
            ylim(ax,[min(y) max(y)]);


            %% Wavy-wall solid region

            fill_patch_wavywall( ...
                x_wall, ...
                y1_wall, ...
                y2_wall, ...
                solid_color);


            %% Wall boundaries

            plot( ...
                x_wall, ...
                y1_wall, ...
                solid_boundary_color, ...
                'LineWidth',solid_boundary_linewidth);

            plot( ...
                x_wall, ...
                y2_wall, ...
                solid_boundary_color, ...
                'LineWidth',solid_boundary_linewidth);

            %% Domain limits
            xlim(ax,[0 Lx]);
            ylim(ax,[min(y) max(y)]);


            %% Ticks

            ax.XTick = linspace( ...
                0, ...
                Lx, ...
                4);

            ax.YTick = linspace( ...
                min(y), ...
                max(y), ...
                5);


            %% Tick labels

            ax.XTickLabel = strip_zeros( ...
                compose('%.1f',ax.XTick));

            ax.YTickLabel = strip_zeros( ...
                compose('%.1f',ax.YTick));


            %% Fonts

            set( ...
                ax, ...
                'FontSize',axis_tick_fontsize);

            set( ...
                cb, ...
                'FontSize',colorbar_tick_fontsize);


            xlabel( ...
                ax, ...
                'x', ...
                'FontSize',axis_label_fontsize);

            ylabel( ...
                ax, ...
                'y', ...
                'FontSize',axis_label_fontsize);
            %% ------------------------------------------------------------
            % EXACT SAME AXES CENTERING AS RESOLVENT POSTPROCESSOR
            % ------------------------------------------------------------
%             set( ...
%                 ax, ...
%                 'Position',[0.22 0.20 0.46 0.54]);
            set(ax,'Position',[0.22 0.28 0.46 0.50]);

            %% Aspect ratio

            daspect([1 1 1]);

            set(fig,'PaperPositionMode','auto');

            %% ========================= SAVE CSV ==========================
            csv_filename = sprintf( ...
                '%s_forcing_c%02d_kz%g.csv', ...
                forc_names{k}, ...
                i, ...
                kz);

            writematrix( ...
                field_plot, ...
                fullfile(snapdir,csv_filename));


            %% ========================= SAVE PNG ==========================

            png_filename = sprintf( ...
                '%s_forcing_c%02d_kz%g.png', ...
                forc_names{k}, ...
                i, ...
                kz);

            exportgraphics( ...
                fig, ...
                fullfile(snapdir,png_filename), ...
                'Resolution',output_resolution, ...
                'BackgroundColor','white');
            close(fig);

        end

    end

end

%% ========================= FINISHED =====================================

fprintf('\nPostprocessing complete.\n');
fprintf('Outputs saved in:\n%s\n',snapdir);

%% ========================================================================
%  LOCAL FUNCTION
%  Same wavy-wall filling logic as postprocess_resolvent_wavywall
%% ========================================================================

function fill_patch_wavywall( ...
    x_wall, ...
    y1_wall, ...
    y2_wall, ...
    solid_color)

ax = gca;

hold(ax,'on');

ymin = ax.YLim(1);
ymax = ax.YLim(2);

xcol = x_wall(:);

y1 = y1_wall(:);
y2 = y2_wall(:);


%% Bottom solid region

xb = [xcol; flipud(xcol)];

yb = [ ...
    y1; ...
    ymin*ones(size(xcol))];

patch( ...
    xb, ...
    yb, ...
    solid_color, ...
    'EdgeColor','none', ...
    'FaceAlpha',1.0, ...
    'Parent',ax);


%% Top solid region

xt = [xcol; flipud(xcol)];

yt = [ ...
    ymax*ones(size(xcol)); ...
    flipud(y2)];

patch( ...
    xt, ...
    yt, ...
    solid_color, ...
    'EdgeColor','none', ...
    'FaceAlpha',1.0, ...
    'Parent',ax);

end


%% ========================================================================
%  LOCAL FUNCTION
%  Remove trailing zeros from tick labels
%% ========================================================================
function out = strip_zeros(lbls)

out = cellstr(lbls);

for k = 1:numel(out)

    s = out{k};

    s = regexprep( ...
        s, ...
        '(\.\d*?)0+$', ...
        '$1');

    s = regexprep( ...
        s, ...
        '\.$', ...
        '');

    out{k} = s;

end

end