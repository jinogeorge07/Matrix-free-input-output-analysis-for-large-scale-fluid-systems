%% postprocess_resolvent.m
% Load resolvent_results.mat and produce plots / saved data for postprocessing

clearvars; close all; clc;

% %% ========================= USER INPUTS ==================================
folderpath = './'
cd(folderpath)

% Results file produced by the resolvent run
%results_file = 'resolvent_results_matrixfree_NxNy116.mat';

% Output folder
out_dir = 'output_figure_3a,b';
plot_log_scale = true;                     % use semilogy for sigma plot
save_mat_processed = false;                 % save processed arrays
%colormap_choice = turbo(256);
colormap_choice = bluewhitered(256);

% --- font size controls for the 2D contour plots ---
big_fontsize = 60;
axis_tick_fontsize = 60;
axis_label_fontsize = 64;
colorbar_tick_fontsize = 60;

% --- create output dir ---
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

% --- grid sizes and indexing ---
Nx = 300;
Ny = 300;
N = Nx*Ny;
Lx = 2*pi;
Ly = 2;
Re = 358;

S1 = load('real_u_forc.mat');
real_x_forc = S1.real_u_forc;

S2 = load('real_v_forc.mat');
real_y_forc = S2.real_v_forc;

S3 = load('real_w_forc.mat');
real_z_forc = S3.real_w_forc;

load('real_u_resp.mat');
load('real_v_resp.mat');
load('real_w_resp.mat');
load('x.mat');
load('y.mat');

% --- 2D contour plots of response components (real part, symmetric colorbar) ---
comp_names = {'u','v','w'};
resp_cells = {real_u_resp, real_v_resp, real_w_resp};
for k = 1:1
    fig = figure('Visible','off','Units','pixels','Position',[100 100 1900 1300]);

    %imagesc(x, y, resp_cells{k}); axis xy; cb = colorbar;
    contourf(x, y, resp_cells{k},40); axis xy; cb = colorbar;
    colormap(colormap_choice);
    % Symmetric color limits about 0, sized to this panel's own max
    % magnitude. Guard against an all-zero field (clim([0 0]) errors).
    clim_val = max(abs(resp_cells{k}(:)));
    if clim_val == 0
        clim_val = eps;
    end
    ax = gca;
    clim(ax,[-clim_val,clim_val]);

    cb = colorbar;
    cb.Ticks = linspace(-clim_val,clim_val,7);
    cb.TickLabels = strip_zeros(compose('%.4f',cb.Ticks));

    colormap(ax,bluewhitered(256));

    % Domain limits
    xlim(ax,[0 Lx]);
    ylim(ax,[min(y) max(y)]);

    % Same tick layout as wavy-wall plots
    ax.XTick = linspace(0,Lx,4);
    ax.YTick = linspace(min(y),max(y),5);

    ax.XTickLabel = strip_zeros(compose('%.1f',ax.XTick));
    ax.YTickLabel = strip_zeros(compose('%.1f',ax.YTick));

    set(ax,'FontSize',axis_tick_fontsize);
    set(cb,'FontSize',colorbar_tick_fontsize);

    xlabel(ax,'x','FontSize',axis_label_fontsize);
    ylabel(ax,'y','FontSize',axis_label_fontsize);

    set(ax,'Position',[0.22 0.20 0.46 0.58]);

    exportgraphics(fig, fullfile(out_dir, sprintf('resp_%s_2D_mf.png', comp_names{k})), ...
        'Resolution', 300);
    saveas(fig, fullfile(out_dir, sprintf('resp_%s_2D_mf.fig', comp_names{k})));
    close(fig);
end

% --- 2D contour plots of forcing components (real part, symmetric colorbar) ---
comp_names = {'x','y','z'};
forc_cells = {real_x_forc, real_y_forc, real_z_forc};
for k = 1:1
    fig = figure('Visible','off','Units','pixels','Position',[100 100 1900 1300]);
    %     imagesc(x, y, forc_cells{k}); axis xy; cb = colorbar;
    contourf(x, y, forc_cells{k},40); axis xy; cb = colorbar;

    colormap(colormap_choice);
    clim_val = max(abs(forc_cells{k}(:)));
    if clim_val == 0
        clim_val = eps;
    end
       ax = gca;
    % --- Same ordering/style as the response loop above. ---
    clim(ax,[-clim_val, clim_val]);

    cb = colorbar;
    cb.Ticks = linspace(-clim_val,clim_val,7);
    cb.TickLabels = strip_zeros(compose('%.4f',cb.Ticks));

    colormap(ax,bluewhitered(256));

    % Domain limits
    xlim(ax,[0 Lx]);
    ylim(ax,[min(y) max(y)]);

    % Same tick layout as wavy-wall plots
    ax.XTick = linspace(0,Lx,4);      % 4 x-ticks including 0 and Lx
    ax.YTick = linspace(min(y),max(y),5);    % 5 y-ticks including top/bottom

    ax.XTickLabel = strip_zeros(compose('%.1f',ax.XTick));
    ax.YTickLabel = strip_zeros(compose('%.1f',ax.YTick));

    set(ax,'FontSize',axis_tick_fontsize);
    set(cb,'FontSize',colorbar_tick_fontsize);

    xlabel(ax,'x','FontSize',axis_label_fontsize);
    ylabel(ax,'y','FontSize',axis_label_fontsize);
    set(ax,'Position',[0.22 0.20 0.46 0.58]);

    exportgraphics(fig, fullfile(out_dir, sprintf('forc_%s_2D_mf.png', comp_names{k})), ...
        'Resolution', 300);
    saveas(fig, fullfile(out_dir, sprintf('forc_%s_2D_mf.fig', comp_names{k})));
    close(fig);
end

% --- optionally write CSVs for quick inspection in Python/R ---
writematrix(real_u_resp, fullfile(out_dir,'real_u_resp.csv'));
writematrix(real_x_forc, fullfile(out_dir,'real_x_forc.csv'));

% --- final message ---
fprintf('Postprocessing complete. Outputs saved in "%s"\n', out_dir);

%% ========================= LOCAL FUNCTION ===============================
% Copied from the input-output postprocessing script so tick labels
% look identical between the two: strips trailing zeros so "0.100"
% becomes "0.1" while multi-digit values keep their needed precision.
function out = strip_zeros(lbls)
out = cellstr(lbls);
for k = 1:numel(out)
    s = out{k};
    s = regexprep(s,'(\.\d*?)0+$','$1');  % drop trailing zeros
    s = regexprep(s,'\.$','');            % drop trailing dot
    out{k} = s;
end
end