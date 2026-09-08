%% postprocess_resolvent.m
% Load resolvent_results.mat and produce plots / saved data for postprocessing

clearvars; close all; clc;

%folderpath = 'E:\input-otuput matrix free results\matlab_26704037';
%folderpath = 'E:\input-otuput matrix free results\compare results\matlab_26862155_compare_matrix';
%folderpath = 'E:\input-otuput matrix free results\matlab_26983385_300x300_matrix_free';
folderpath = 'E:\input-otuput matrix free results\symmetric wavywall\matlab_43905321_400x400_matrixfree_kz1.0';

%E:\input-otuput matrix free results\matlab_26983385_300x300_matrix_free
% C:\Users\jinog\Documents\MATLAB
cd(folderpath)


% --- PARAMETERS you can edit ---
results_file = 'resolvent_results_matrixfree_NxNy400.mat';   % file produced by your run
%results_file = 'resolvent_results_matrixfree_NxNy300';   % file produced by your run

out_dir = 'postproc_output_mf_400x400';               % where figures and data will be saved
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

% --- load results ---
if ~isfile(results_file)
    error('Results file not found: %s', results_file);
end
S = load(results_file);

% % Required variables check
required = {'U_mf','S_mf','V_mf','params','sigma_mf'};
for k = 1:numel(required)
    if ~isfield(S, required{k})
        error('Missing variable "%s" in %s', required{k}, results_file);
    end
end


U_mf = S.U_mf;
S_mf = S.S_mf;
V_mf = S.V_mf;
params = S.params;
sigma_mf = S.sigma_mf;

% --- grid sizes and indexing ---
Nx = params.Nx;
Ny = params.Ny;
N = Nx*Ny;
if size(U_mf,1) ~= 3*N
    error('Unexpected U_mf size: expected 3*N rows (3*%d), got %d', N, size(U_mf,1));
end

% --- singular values ---
singvals = diag(S_mf);
nsv = numel(singvals);

% Plot singular values
fig = figure('Visible','off');
if plot_log_scale
    semilogy(1:nsv, singvals, 'o-', 'LineWidth', 1.5);
else
    plot(1:nsv, singvals, 'o-', 'LineWidth', 1.5);
end
xlabel('singular value index'); ylabel('\sigma'); grid on;
title('Singular values (matrix-free)');
saveas(fig, fullfile(out_dir,'singular_values.png'));
saveas(fig, fullfile(out_dir,'singular_values.fig'));
close(fig);

% --- reshape left singular vector (response) and right singular vector (forcing) ---
% ordering in code: [u; v; w] each of length N (Ny*Nx)
Uvec = U_mf(:,1).*params.w_all.^(-1/2);
Vvec = V_mf(:,1).*params.w_all.^(-1/2);


u_resp = reshape(Uvec(1:N), Ny, Nx);
v_resp = reshape(Uvec(N+1:2*N), Ny, Nx);
w_resp = reshape(Uvec(2*N+1:3*N), Ny, Nx);

u_forc = reshape(Vvec(1:N), Ny, Nx);
v_forc = reshape(Vvec(N+1:2*N), Ny, Nx);
w_forc = reshape(Vvec(2*N+1:3*N), Ny, Nx);

% compute absolute magnitudes (kept for the wall-normal profile/slice
% plots below, and for the energy-norm summary -- magnitude is the
% right quantity there, it's only the 2D colorbar plots that needed
% signed data instead).
abs_u_resp = abs(u_resp);
abs_v_resp = abs(v_resp);
abs_w_resp = abs(w_resp);

abs_u_forc = abs(u_forc);
abs_v_forc = abs(v_forc);
abs_w_forc = abs(w_forc);


real_u_resp = real(u_resp);
real_v_resp = real(v_resp);
real_w_resp = real(w_resp);

real_u_forc = real(u_forc);
real_v_forc = real(v_forc);
real_w_forc = real(w_forc);

% --- coordinate vectors (y from cheb points, x uniform Fourier) ---
y = params.cheb_y;            % cheb points in [-1,1] (size Ny)
% build x vector consistent with fourdif scaling (0..Lx)
%x = linspace(0, params.Lx, Nx+1); x = x(1:end-1);  % Nx points periodic
x = linspace(0, params.Lx, Nx+1); x = x(1:end-1);  % Nx points periodic

% --- 2D contour plots of response components (real part, symmetric colorbar) ---
comp_names = {'u','v','w'};
resp_cells = {real_u_resp, real_v_resp, real_w_resp};
for k = 1:3
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

    exportgraphics(fig, fullfile(out_dir, sprintf('resp_%s_2D_mf.png', comp_names{k})), ...
        'Resolution', 300);
    saveas(fig, fullfile(out_dir, sprintf('resp_%s_2D_mf.fig', comp_names{k})));
    close(fig);
end

% --- 2D contour plots of forcing components (real part, symmetric colorbar) ---
forc_cells = {real_u_forc, real_v_forc, real_w_forc};
for k = 1:3
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

    exportgraphics(fig, fullfile(out_dir, sprintf('forc_%s_2D_mf.png', comp_names{k})), ...
        'Resolution', 300);
    saveas(fig, fullfile(out_dir, sprintf('forc_%s_2D_mf.fig', comp_names{k})));
    close(fig);
end

% --- x-averaged wall-normal profiles (|.| averaged over x) ---
% Unchanged: magnitude is the right quantity for a profile plot (no
% colorbar involved here, so the symmetric-about-0 fix doesn't apply).
ux_profile = mean(abs_u_resp,2);
vx_profile = mean(abs_v_resp,2);
wx_profile = mean(abs_w_resp,2);

fig = figure('Visible','off');
plot(ux_profile, y, '-o', 'LineWidth', 1.2); hold on;
plot(vx_profile, y, '-s', 'LineWidth', 1.2);
plot(wx_profile, y, '-^', 'LineWidth', 1.2);
xlabel('x-averaged |component|'); ylabel('y'); grid on;
legend('|u|','|v|','|w|','Location','best');
title('Wall-normal profiles (x-averaged) of response');
saveas(fig, fullfile(out_dir,'resp_profiles_xavg.png'));
saveas(fig, fullfile(out_dir,'resp_profiles_xavg.fig'));
close(fig);

% --- x-slice at center (x index) ---
ix_center = ceil(Nx/2);
fig = figure('Visible','off');
plot(y, abs_u_resp(:,ix_center), '-o', 'LineWidth', 1.2); hold on;
plot(y, abs_v_resp(:,ix_center), '-s', 'LineWidth', 1.2);
plot(y, abs_w_resp(:,ix_center), '-^', 'LineWidth', 1.2);
xlabel('y'); ylabel('|component| at x_{center}'); grid on;
legend('|u|','|v|','|w|','Location','best');
title(sprintf('Wall-normal slice at x index %d', ix_center));
saveas(fig, fullfile(out_dir,'resp_slice_xcenter.png'));
saveas(fig, fullfile(out_dir,'resp_slice_xcenter.fig'));
close(fig);

% --- energy norm of response and forcing (global) ---
% Unchanged: these are norms, magnitude (abs) is the correct quantity.
E_resp = sum(abs(Uvec).^2);
E_forc = sum(abs(Vvec).^2);
fid = fopen(fullfile(out_dir,'energy_summary.txt'),'w');
fprintf(fid,'Largest singular value (sigma): %g\n', singvals(1));
fprintf(fid,'Response energy (||U||^2): %g\n', E_resp);
fprintf(fid,'Forcing energy (||V||^2): %g\n', E_forc);
fclose(fid);

% --- save processed arrays for later postprocessing ---
if save_mat_processed
    save(fullfile(out_dir,'processed_fields.mat'), ...
        'u_resp','v_resp','w_resp', ...
        'u_forc','v_forc','w_forc', ...
        'abs_u_resp','abs_v_resp','abs_w_resp', ...
        'abs_u_forc','abs_v_forc','abs_w_forc', ...
        'real_u_resp','real_v_resp','real_w_resp', ...
        'real_u_forc','real_v_forc','real_w_forc', ...
        'x','y','singvals','sigma_mf','params','-v7.3');
end

% --- optionally write CSVs for quick inspection in Python/R ---
% Kept the original abs(...) exports so nothing existing breaks, and
% added real(...) exports since that's what the 2D plots now show.
writematrix(abs_u_resp, fullfile(out_dir,'abs_u_resp.csv'));
writematrix(abs_v_resp, fullfile(out_dir,'abs_v_resp.csv'));
writematrix(abs_w_resp, fullfile(out_dir,'abs_w_resp.csv'));

writematrix(real_u_resp, fullfile(out_dir,'real_u_resp.csv'));
writematrix(real_v_resp, fullfile(out_dir,'real_v_resp.csv'));
writematrix(real_w_resp, fullfile(out_dir,'real_w_resp.csv'));

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