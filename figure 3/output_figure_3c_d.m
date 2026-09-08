%% postprocess_resolvent.m
% Load resolvent_results.mat and produce plots / saved data for postprocessing

clearvars; close all; clc;
% %% ========================= USER INPUTS ==================================
folderpath = './'
cd(folderpath)

% --- PARAMETERS you can edit ---
results_file = 'resolvent_results_full_120.mat';   % file produced by your run
out_dir = 'output_figure_3c_d';               % where figures and data will be saved

plot_log_scale = true;                     % use semilogy for sigma plot
save_mat_processed = true;                 % save processed arrays
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

% --- load results ---
if ~isfile(results_file)
    error('Results file not found: %s', results_file);
end
S = load(results_file);

% Required variables check
required = {'U_full_i','S_full_i','V_full_i','params_i','sigma_full'};
for k = 1:numel(required)
    if ~isfield(S, required{k})
        error('Missing variable "%s" in %s', required{k}, results_file);
    end
end

U_full = S.U_full_i;
S_full = S.S_full_i;
V_full = S.V_full_i;
params = S.params_i;
sigma_full = S.sigma_full;

% --- grid sizes and indexing ---
Nx = params.Nx;
Ny = params.Ny;
N = Nx*Ny;
if size(U_full,1) ~= 3*N
    error('Unexpected U_full size: expected 3*N rows (3*%d), got %d', N, size(U_full,1));
end

Uvec = U_full(:,1).*params.w_all.^(-1/2);
Vvec =  V_full(:,1).*params.w_all.^(-1/2);

u_resp = reshape(Uvec(1:N), Ny, Nx);
v_resp = reshape(Uvec(N+1:2*N), Ny, Nx);
w_resp = reshape(Uvec(2*N+1:3*N), Ny, Nx);

u_forc = reshape(Vvec(1:N), Ny, Nx);
v_forc = reshape(Vvec(N+1:2*N), Ny, Nx);
w_forc = reshape(Vvec(2*N+1:3*N), Ny, Nx);

abs_u_resp = abs(u_resp);
abs_v_resp = abs(v_resp);
abs_w_resp = abs(w_resp);

abs_u_forc = abs(u_forc);
abs_v_forc = abs(v_forc);
abs_w_forc = abs(w_forc);

% colormap is actually for.
real_u_resp = real(u_resp);
real_v_resp = real(v_resp);
real_w_resp = real(w_resp);

real_u_forc = real(u_forc);
real_v_forc = real(v_forc);
real_w_forc = real(w_forc);

% --- coordinate vectors (y from cheb points, x uniform Fourier) ---
y = params.cheb_y;            % cheb points in [-1,1] (size Ny)
% build x vector consistent with fourdif scaling (0..Lx)
x = linspace(0, params.Lx, Nx+1); x = x(1:end-1);  % Nx points periodic

% --- 2D contour plots of response components (real part, symmetric colorbar) ---
comp_names = {'u','v','w'};
resp_cells = {real_u_resp, real_v_resp, real_w_resp};
for k = 1:3
    fig = figure('Visible','off','Units','pixels','Position',[100 100 1900 1300]);

    %contourf(x, y, resp_cells{k},40); axis xy;
    contourf(x, y, resp_cells{k},40); axis xy; cb = colorbar;

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
    xlim(ax,[0 params.Lx]);
    ylim(ax,[min(y) max(y)]);

    % Same tick layout as wavy-wall plots
    ax.XTick = linspace(0,params.Lx,4);
    ax.YTick = linspace(min(y),max(y),5);

    ax.XTickLabel = strip_zeros(compose('%.1f',ax.XTick));
    ax.YTickLabel = strip_zeros(compose('%.1f',ax.YTick));

    set(ax,'FontSize',axis_tick_fontsize);
    set(cb,'FontSize',colorbar_tick_fontsize);

    xlabel(ax,'x','FontSize',axis_label_fontsize);
    ylabel(ax,'y','FontSize',axis_label_fontsize);

    set(ax,'Position',[0.22 0.20 0.46 0.58]);
    %title(sprintf('Re(%s) (response) — largest singular vector', comp_names{k}));

    % --- FIX: exportgraphics instead of saveas for reliable cropping. ---
    % exportgraphics recomputes the actual visible bounding box at save
    % time (including axis/colorbar labels and tick text), rather than
    % using a fixed crop the way saveas does -- more robust once font
    % sizes get this large.
    exportgraphics(fig, fullfile(out_dir, sprintf('resp_%s_2D.png', comp_names{k})), ...
        'Resolution', 300);
    saveas(fig, fullfile(out_dir, sprintf('resp_%s_2D.fig', comp_names{k})));
    close(fig);
end

% --- 2D contour plots of forcing components (real part, symmetric colorbar) ---
comp_names = {'x','y','z'};
forc_cells = {real_u_forc, real_v_forc, real_w_forc};
for k = 1:3
    fig = figure('Visible','off','Units','pixels','Position',[100 100 1900 1300]);

    %contourf(x, y, forc_cells{k},40); axis xy;
    contourf(x, y, forc_cells{k},40); axis xy; cb = colorbar;

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
    xlim(ax,[0 params.Lx]);
    ylim(ax,[min(y) max(y)]);

    % Same tick layout as wavy-wall plots
    ax.XTick = linspace(0,params.Lx,4);      % 4 x-ticks including 0 and Lx
    ax.YTick = linspace(min(y),max(y),5);    % 5 y-ticks including top/bottom

    ax.XTickLabel = strip_zeros(compose('%.1f',ax.XTick));
    ax.YTickLabel = strip_zeros(compose('%.1f',ax.YTick));

    set(ax,'FontSize',axis_tick_fontsize);
    set(cb,'FontSize',colorbar_tick_fontsize);

    xlabel(ax,'x','FontSize',axis_label_fontsize);
    ylabel(ax,'y','FontSize',axis_label_fontsize);

    set(ax,'Position',[0.22 0.20 0.46 0.58]);

    exportgraphics(fig, fullfile(out_dir, sprintf('forc_%s_2D.png', comp_names{k})), ...
        'Resolution', 300);
    saveas(fig, fullfile(out_dir, sprintf('forc_%s_2D.fig', comp_names{k})));
    close(fig);
end

% --- optionally write CSVs---

writematrix(real_u_resp, fullfile(out_dir,'real_u_resp.csv'));
writematrix(real_u_forc, fullfile(out_dir,'real_x_forc.csv'));

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