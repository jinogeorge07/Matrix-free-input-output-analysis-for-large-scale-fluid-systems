%% postprocess_resolvent_wavywall.m
clearvars; close all; clc;

% %% ========================= USER INPUTS ==================================
folderpath = './'
cd(folderpath)

% Results file produced by the resolvent run
%results_file = 'resolvent_results_matrixfree_NxNy116.mat';

% Output folder
out_dir = 'output_figure_6a,b';

% Plot controls
plot_log_scale = true;
colormap_choice = bluewhitered(256);

% Font size for the 2D plots
big_fontsize = 52;
axis_tick_fontsize = 60;
axis_label_fontsize = 64;
colorbar_tick_fontsize = 60;

% Solid-region display
solid_color = [0.30 0.30 0.30];
solid_boundary_color = 'k';
solid_boundary_linewidth = 2.0;

epsilon = 0.12;

% load('U_mf.mat');
% load('V_mf.mat');
% load('sigma_mf.mat');
% load('params.mat');

Nx = 116;
Ny = 156;
N = Nx*Ny;
Re = 500;

%% ========================= LOAD GEOMETRY ================================
load('x_solid.mat');
load('y_solid.mat');

fname = sprintf( ...
    'mask_smooth_wavywall_Re%d_%dx%d_epsilon_%.2f.mat', ...
    Re,Ny,Nx,epsilon);

load(fullfile(fname));

solid_mask = (mask_smooth == Re);

%% ========================= OUTPUT FOLDER ================================
%out_dir = fullfile(out_dir);

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

load('real_u_forc.mat');
load('real_v_forc.mat');
load('real_w_forc.mat');
load('real_u_resp.mat');
load('real_v_resp.mat');
load('real_w_resp.mat');
load('x.mat');
load('y.mat');

%% ========================= COORDINATES ==================================
% y = params.cheb_y;
% x = params.x;
%x = x(:).';
%y = y(:);

Lx = 2.4;
y0 = 1;
A1 = epsilon;
A2 = epsilon;

%% ------------------------------------------------------------------------
% Wavy-wall curves on the same x-grid as the contour plots
x_wall  = x(:).';
y1_wall = -y0 - A1*sin(2*pi/Lx*x_wall);
y2_wall =  y0 + A2*sin(2*pi/Lx*x_wall);
%% ------------------------------------------------------------------------
y_axis_limit = 1.20;

%% ========================= RESPONSE PLOTS ===============================
comp_names = {'u','v','w'};
resp_cells = {real_u_resp,real_v_resp,real_w_resp};

for k = 1:1

    field_plot = resp_cells{k};

%     fig = figure( ...
%         'Visible','off', ...
%         'Units','pixels', ...
%         'Position',[100 100 1900 1300]); % 100 100 1350 1000

fig = figure( ...
        'Visible','off', ...
        'Units','pixels', ...
        'Position',[100 100 1350 1000]); % 

    contourf(x,y,field_plot,40);
    axis xy;
    hold on;

    ylim(gca,[-y_axis_limit,y_axis_limit]);

    %clim_val = max(abs(field_plot(:)));
    clim_val = 0.173;
    ax = gca;

    clim(ax,[-clim_val,clim_val]);

    cb = colorbar;
    cb.Ticks = linspace(-clim_val,clim_val,7);
    cb.TickLabels = strip_zeros(compose('%.3f',cb.Ticks));

    colormap(ax,bluewhitered(256));

    % Cover the regions below and above the wavy walls
    fill_patch_wavywall(x_wall,y1_wall,y2_wall,solid_color);

    % Plot the wall boundaries
    plot(x_wall,y1_wall,solid_boundary_color, ...
        'LineWidth',solid_boundary_linewidth);

    plot(x_wall,y2_wall,solid_boundary_color, ...
        'LineWidth',solid_boundary_linewidth);

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

    %set(gca,'Position',[0.22 0.20 0.46 0.58]);
    set(gca,'Position',[0.22 0.28 0.46 0.50]);

    daspect([1 1 1]);
    set(fig,'PaperPositionMode','auto');

    exportgraphics(fig, ...
    fullfile(out_dir,sprintf('resp_%s_2D_mf_ww.png',comp_names{k})), ...
    'Resolution',300);

%     print(fig, ...
%         fullfile(out_dir,sprintf('resp_%s_2D_mf',comp_names{k})), ...
%         '-dpng','-r300');

%     saveas(fig, ...
%         fullfile(out_dir,sprintf('resp_%s_2D_mf_ww.fig',comp_names{k})));
    close(fig);
end

%% ========================= FORCING PLOTS ================================
forc_cells = {real_u_forc,real_v_forc,real_w_forc};

for k = 1:1

    field_plot = forc_cells{k};

%     fig = figure( ...
%         'Visible','off', ...
%         'Units','pixels', ...
%         'Position',[100 100 1900 1300]);

fig = figure( ...
        'Visible','off', ...
        'Units','pixels', ...
        'Position',[100 100 1350 1000]);

    contourf(x,y,field_plot,40);
    axis xy;
    hold on;

    % Set ylim BEFORE fill_patch_wavywall runs -- same reason as the
    % response loop above.
    ylim(gca,[-y_axis_limit,y_axis_limit]);

    %clim_val = max(abs(field_plot(:)));
    clim_val = 0.151;
    ax = gca;

    % --- Same ordering/style as the response loop above. ---
    clim(ax,[-clim_val,clim_val]);

    cb = colorbar;
    cb.Ticks = linspace(-clim_val,clim_val,7);
    cb.TickLabels = strip_zeros(compose('%.3f',cb.Ticks));

    colormap(ax,bluewhitered(256));

    % Cover the regions below and above the wavy walls
    fill_patch_wavywall(x_wall,y1_wall,y2_wall,solid_color);

    % Plot the wall boundaries
    plot(x_wall,y1_wall,solid_boundary_color, ...
        'LineWidth',solid_boundary_linewidth);

    plot(x_wall,y2_wall,solid_boundary_color, ...
        'LineWidth',solid_boundary_linewidth);

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

    %set(gca,'Position',[0.22 0.20 0.46 0.58]);
     set(gca,'Position',[0.22 0.28 0.46 0.50]);

    %xlabel('x');
    %ylabel('y');

    daspect([1 1 1]);
    set(fig,'PaperPositionMode','auto');

%     print(fig, ...
%         fullfile(out_dir,sprintf('forc_%s_2D_mf',comp_names{k})), ...
%         '-dpng','-r300');
    
      exportgraphics(fig, ...
          fullfile(out_dir,sprintf('forc_%s_2D_mf_ww.png',comp_names{k})), ...
            'Resolution',300);

%     saveas(fig, ...
%         fullfile(out_dir,sprintf('forc_%s_2D_mf_ww.fig',comp_names{k})));

    close(fig);
end


%% ========================= CSV OUTPUT ===================================
writematrix(real_u_resp,fullfile(out_dir,'real_u_resp.csv'));
writematrix(real_u_forc,fullfile(out_dir,'real_u_forc.csv'));

% writematrix(real_v_resp,fullfile(out_dir,'real_v_resp.csv'));
% writematrix(real_w_resp,fullfile(out_dir,'real_w_resp.csv'));

fprintf('Postprocessing complete. Outputs saved in "%s"\n',out_dir);

%% ========================= LOCAL FUNCTION ===============================
function fill_patch_wavywall(x_wall,y1_wall,y2_wall,solid_color)

ax = gca;
hold(ax,'on');

ymin = ax.YLim(1);
ymax = ax.YLim(2);

xcol = x_wall(:);
y1 = y1_wall(:);
y2 = y2_wall(:);

% Region below the bottom wall
xb = [xcol; flipud(xcol)];
yb = [y1; ymin*ones(size(xcol))];

patch(xb,yb,solid_color, ...
    'EdgeColor','none', ...
    'FaceAlpha',1.0, ...
    'Parent',ax);

% Region above the top wall
xt = [xcol; flipud(xcol)];
yt = [ymax*ones(size(xcol)); flipud(y2)];

patch(xt,yt,solid_color, ...
    'EdgeColor','none', ...
    'FaceAlpha',1.0, ...
    'Parent',ax);

end

%% ========================= LOCAL FUNCTION ===============================
% Copied from the input-output postprocessing script so tick labels
% look identical across all three postprocessing scripts: strips
% trailing zeros so "0.100" becomes "0.1" while multi-digit values
% keep their needed precision.
function out = strip_zeros(lbls)
out = cellstr(lbls);
for k = 1:numel(out)
    s = out{k};
    s = regexprep(s,'(\.\d*?)0+$','$1');  % drop trailing zeros
    s = regexprep(s,'\.$','');            % drop trailing dot
    out{k} = s;
end
end