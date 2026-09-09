clearvars;
close all;
clc;

% %% ========================= USER INPUTS ==================================
folderpath = './'
cd(folderpath)

load data_x.mat
load data_y.mat
load u_mean_zt.mat
load v_mean_zt.mat
load dUdx_mean_zt.mat
load dUdy_mean_zt.mat
load dVdx_mean_zt.mat
load dVdy_mean_zt.mat

Ny=30;
Nx=30;
epsilon = 0.12;
h = 1; y0 = h;
A1 = epsilon; A2 = epsilon;
Lx = 2.4
Ly1 = 1.0 + 1.0999*epsilon; % Half-height
Re = 500;

%mask_smooth_wavywall_Re500_156x116_epsilon_0.12
fname = sprintf('mask_smooth_wavywall_Re%d_%dx%d_epsilon_%.2f.mat', Re, Ny, Nx, epsilon);
load(fname);

x_old = data_x; y_old = data_y;

%% -------------fine mesh from Dedalus-------------------------------------
[X_old, Y_old] = meshgrid(x_old, y_old);
%% ------------------------------------------------------------------------

% (1) Robust interpolation setup: Dedalus fields often [Nx x Ny] → transpose
U_f  = u_mean_zt;  V_f  = v_mean_zt;     % [Ny_f x Nx_f]
Ux_f = dUdx_mean_zt; Uy_f = dUdy_mean_zt;
Vx_f = dVdx_mean_zt; Vy_f = dVdy_mean_zt;
K_f  = mask_smooth;                        % [Ny_f x Nx_f]

%% Spectral operators
[x,Dxx] = fourdif(Nx,2);
[~,Dx]  = fourdif(Nx,1);
[cheb_y,DM] = chebdif(Ny,2);
D1 = DM(:,:,1);
D2 = DM(:,:,2);

% Scaling
% D1  = -D1*(2/(2+2*epsilon)); % negative bc we flipped y so derivative also needs to be flipped
% D2  = D2*(2/(2+2*epsilon))^2;

Dy_mat  = -D1*(2/(2*Ly1)); % negative bc we flipped y so derivative also needs to be flipped
Dyy_mat  = D2*(2/(2*Ly1))^2;
Dx_mat  = Dx*(2*pi/Lx);
Dxx_mat = Dxx*(2*pi/Lx)^2;

Ix = speye(Nx); Iy = speye(Ny);
I  = kron(Ix,Iy);  % sparse because Ix, Iy are sparse

[~,w] = clencurt(Ny-1);
%w = kron(Ix, diag(w));     % sparse
% x = x * (Lx / (2*pi));
% y = cheb_y/2*(2+2*1.1*epsilon);
% y = flip(y); % chebyshev is in descending order so needs to be flipped
% w = w*((2+2*epsilon)/2); % normalize weighted matrix

x = x * (Lx / (2*pi));

y = cheb_y/2*(2*Ly1);
y = flip(y); % chebyshev is in descending order so needs to be flipped
w = w*((2*Ly1)/2); % normalize weighted matrix

w_sqrt = sqrt(w);
w_sqrt = spdiags(w_sqrt',0,Ny,Ny);
w_sqrt = kron(Ix, w_sqrt);

w_sqrt_inv = 1./(sqrt(w)); % elementwise sqrt
w_sqrt_inv = spdiags(w_sqrt_inv',0,Ny,Ny);
w_sqrt_inv = kron(Ix, w_sqrt_inv);

%% Coarse grid
[X_new, Y_new] = meshgrid(x, y);

%% (1) Interpolate (vector form) to coarse grid
U_xy_coarse   = interp2(X_old, Y_old, U_f,  X_new, Y_new, 'linear');
V_xy_coarse   = interp2(X_old, Y_old, V_f,  X_new, Y_new, 'linear');
dUdx_xy_coarse= interp2(X_old, Y_old, Ux_f, X_new, Y_new, 'linear');
dUdy_xy_coarse= interp2(X_old, Y_old, Uy_f, X_new, Y_new, 'linear');
dVdx_xy_coarse= interp2(X_old, Y_old, Vx_f, X_new, Y_new, 'linear');
dVdy_xy_coarse= interp2(X_old, Y_old, Vy_f, X_new, Y_new, 'linear');

%% Flattened fields (keep your names)
N = Ny*Nx;
c_number = 12;

% params.kx=0;
params.kz=1;
params.omega=1;
params.Ny=Ny;
params.Nx=Nx;
params.N=params.Nx*params.Ny;
params.Re=Re;
params.Ly=2*Ly1;
params.Lx=Lx;
params.gmres_tol=1e-7;
params.gmres_restart=250;
params.gmres_maxit=params.Nx*params.Ny;
params.Dy_mat = Dy_mat;
params.Dyy_mat = Dyy_mat;
params.Dx_mat = Dx_mat;
params.Dxx_mat = Dxx_mat;
params.cheb_y = y;
params.x = x;

params.U_mat=U_xy_coarse;
params.V_mat=V_xy_coarse;
params.dUdx_mat=dUdx_xy_coarse;
params.dUdy_mat=dUdy_xy_coarse;
params.dVdx_mat=dVdx_xy_coarse;
params.dVdy_mat=dVdy_xy_coarse;
params.K_inv_mat=mask_smooth;

params.U     = U_xy_coarse(:);
params.V     = V_xy_coarse(:);
params.dUdx  = dUdx_xy_coarse(:);
params.dUdy  = dUdy_xy_coarse(:);
params.dVdx  = dVdx_xy_coarse(:);
params.dVdy  = dVdy_xy_coarse(:);
params.K_inv = mask_smooth(:);

U_xave     = mean(U_xy_coarse,2);
V_xave     = mean(V_xy_coarse,2);
dUdx_xave  = mean(dUdx_xy_coarse,2);
dUdy_xave  = mean(dUdy_xy_coarse,2);
dVdx_xave  = mean(dVdx_xy_coarse,2);
dVdy_xave  = mean(dVdy_xy_coarse,2);
Kinv_xave  = mean(mask_smooth,2);

params.U_xave_diag     = spdiags(U_xave,0,Ny,Ny);
params.V_xave_diag     = spdiags(V_xave,0,Ny,Ny);
params.dUdx_xave_diag  = spdiags(dUdx_xave,0,Ny,Ny);
params.dUdy_xave_diag  = spdiags(dUdy_xave,0,Ny,Ny);
params.dVdx_xave_diag  = spdiags(dVdx_xave,0,Ny,Ny);
params.dVdy_xave_diag  = spdiags(dVdy_xave,0,Ny,Ny);
params.K_inv_xave_diag = spdiags(Kinv_xave,0,Ny,Ny);

%% if compare_full_matrix is on, then you have to turn on matrix_free_only as well
% matrix_free = false;
% compare_full_matrix= false;
% run_grid_benchmark = false;
% full_matrix_only = true;

%% ================= RUN OPTIONS =================
run_grid_benchmark  = false;   % benchmark multiple grids: full vs matrix-free
matrix_free         = true;   % matrix-free calculation for current grid
full_matrix_only    = false;    % full-matrix calculation only for current grid
compare_full_matrix = false;   % compare full vs matrix-free for current grid

% Preconditioner options:
%   'fft_1d_laplacian'    : FFT in x and 1D Chebyshev Laplacian per mode
%   'fft_1d_resolvent'    : FFT in x and coupled 1D resolvent/Stokes block per mode
params.preconditioner_type='fft_1d_resolvent'; %fft_1d_resolvent

if run_grid_benchmark
    params.epsilon = epsilon;

    params.benchmark_source.x_old = x_old;
    params.benchmark_source.y_old = y_old;

    params.benchmark_source.U_f  = U_f;
    params.benchmark_source.V_f  = V_f;
    params.benchmark_source.Ux_f = Ux_f;
    params.benchmark_source.Uy_f = Uy_f;
    params.benchmark_source.Vx_f = Vx_f;
    params.benchmark_source.Vy_f = Vy_f;
end


params=prepare_params(params);

assert(isequal(size(params.U_mat),[Ny,Nx]));
assert(isequal(size(params.V_mat),[Ny,Nx]));
assert(isequal(size(params.K_inv_mat),[Ny,Nx]));

assert(isequal(size(params.U_xave_diag),[Ny,Ny]));
assert(isequal(size(params.V_xave_diag),[Ny,Ny]));
assert(isequal(size(params.K_inv_xave_diag),[Ny,Ny]));

assert(all(isfinite(params.U_mat),'all'), ...
    'U interpolation contains NaN or Inf values.');
assert(all(isfinite(params.V_mat),'all'), ...
    'V interpolation contains NaN or Inf values.');


benchmark_Nx_list = [12 24 36 48 80 120 140];
benchmark_Ny_list = [12 24 36 48 80 120 140];

nsv = 1;
opts.tol = 1e-3;
opts.maxit = 30;
opts.p = 8;
opts.disp=1;

if matrix_free

    % Shared with H_fun; initialize before the main SVDS run.
    % One row per GMRES call, in execution order (forward and adjoint).
    global gmres_info
    gmres_info = table('Size',[0 7], ...
    'VariableTypes',{'double','cell','double','double','double','double','double'}, ...
    'VariableNames',{'solve_index','direction','outer_iterations', ...
    'inner_iterations','total_inner_iterations','flag','relres'});


    matrix_free_total_cpu_start = cputime;
    matrix_free_setup_cpu_start = cputime;
    params = add_laplacian_preconditioner(params);
    matrix_free_setup_cpu_time_s = cputime-matrix_free_setup_cpu_start;
    H_mf = @(f,tflag) H_fun(f,tflag,params);

    matrix_free_svd_cpu_start = cputime; % command to measure cput time
    [U_mf,S_mf,V_mf] = svds(H_mf,[3*Nx*Ny,3*Nx*Ny],nsv,'largest',opts);
    matrix_free_svd_cpu_time_s = cputime-matrix_free_svd_cpu_start;
    matrix_free_total_cpu_time_s = cputime-matrix_free_total_cpu_start;
    sigma_mf = diag(S_mf);
    disp(table(sigma_mf,'VariableNames',{'sigma_matrix_free'}))
    gmres_average_iterations = mean(gmres_info{:,{'inner_iterations','outer_iterations','total_inner_iterations'}},1); % [mean inner, mean outer]

    disp(table(gmres_average_iterations, ...
    'VariableNames', {'mean_inner_outer_total inner'}));

    matrix_free_memory_MB_a = (workspace_bytes( ...
        params,U_mf,S_mf,V_mf)+matrix_free_workspace_bytes(params))/1024^2;

    save(sprintf('resolvent_results_matrixfree_NxNy%d.mat', Nx), ...
        'U_mf', 'S_mf', 'V_mf', ...
        'sigma_mf', ...
        'params', ...
        'matrix_free_setup_cpu_time_s', ...
        'matrix_free_svd_cpu_time_s', ...
        'matrix_free_total_cpu_time_s', ...
        '-v7.3');

        cpu_time_table = table(matrix_free_memory_MB_a, matrix_free_setup_cpu_time_s, ...
        matrix_free_svd_cpu_time_s,matrix_free_total_cpu_time_s, ...
        'VariableNames',{'matrix_free_memory_MB','matrix_free_preconditioner_setup_cpu_s', ...
        'matrix_free_svds_cpu_s','matrix_free_total_cpu_s', ...
        });
    disp(cpu_time_table)
    writetable(cpu_time_table,'cpu_time_table_matrixfree_1D_resolvent_preconditioner.csv');
end

if compare_full_matrix
    full_matrix_total_cpu_start = cputime;
    full_matrix_build_cpu_start = cputime;
    [H,L]=resolvent_2D_xy(params);
    full_matrix_build_cpu_time_s = cputime-full_matrix_build_cpu_start;

    %% Compare full and matrix-free forward L actions
    Ntest = params.N;

    q_test = randn(4*Ntest,1) + 1i*randn(4*Ntest,1);

    Lq_mf   = L_fun(q_test,'notransp',params);
    Lq_full = L*q_test;

    L_action_error = norm(Lq_mf-Lq_full)/norm(Lq_full);

    fprintf('Forward L-action error = %.3e\n',L_action_error);

    %% Compare full and matrix-free forward H actions
    f_test = randn(3*Ntest,1) + 1i*randn(3*Ntest,1);

    Hf_mf   = H_fun(f_test,'notransp',params);
    Hf_full = H*f_test;

    H_action_error = norm(Hf_mf-Hf_full)/norm(Hf_full);

    fprintf('Forward H-action error = %.3e\n',H_action_error);

    full_matrix_svd_cpu_start = cputime;
    [U_full,S_full,V_full] = svds(H,nsv,'largest');
    full_matrix_svd_cpu_time_s = cputime-full_matrix_svd_cpu_start;
    full_matrix_total_cpu_time_s = cputime-full_matrix_total_cpu_start;
    sigma_full = diag(S_full);
    sigma_rel_err = abs(sigma_mf-sigma_full)./abs(sigma_full);

    [U_mf_aligned,u_err] = align_singular_vectors(U_mf,U_full);
    [V_mf_aligned,v_err] = align_singular_vectors(V_mf,V_full);

    comparison_table = table(sigma_mf,sigma_full,sigma_rel_err,u_err.',v_err.', ...
        'VariableNames',{'sigma_matrix_free','sigma_full','sigma_relative_error', ...
        'left_vector_error','right_vector_error'});
    disp(comparison_table)

    cpu_time_table = table(matrix_free_setup_cpu_time_s, ...
        matrix_free_svd_cpu_time_s,matrix_free_total_cpu_time_s, ...
        full_matrix_build_cpu_time_s,full_matrix_svd_cpu_time_s, ...
        full_matrix_total_cpu_time_s, ...
        'VariableNames',{'matrix_free_preconditioner_setup_cpu_s', ...
        'matrix_free_svds_cpu_s','matrix_free_total_cpu_s', ...
        'full_matrix_build_cpu_s','full_matrix_svds_cpu_s', ...
        'full_matrix_total_cpu_s'});
    disp(cpu_time_table)

    figure;
    tiledlayout(2,2)

    nexttile
    bar(categorical({'matrix-free','full matrix'}),[sigma_mf(1),sigma_full(1)])
    ylabel('$\sigma$','Interpreter','latex')
    title('Largest singular value')

    nexttile
    semilogy(1:nsv,sigma_rel_err,'o-')
    xlabel('singular value index')
    ylabel('relative error')
    title('Singular value error')
    grid on

    nexttile
    plot(abs(U_full(:,1)),'k-','LineWidth',1.5);
    hold on
    plot(abs(U_mf_aligned(:,1)),'r--','LineWidth',1.2)
    legend('full matrix','matrix-free','Location','best')
    xlabel('state index')
    ylabel('$|u_1|$','Interpreter','latex')
    title('Left singular vector')
    grid on

    nexttile
    plot(abs(V_full(:,1)),'k-','LineWidth',1.5)
    hold on
    plot(abs(V_mf_aligned(:,1)),'r--','LineWidth',1.2)
    legend('full matrix','matrix-free','Location','best')
    xlabel('forcing index')
    ylabel('$|v_1|$','Interpreter','latex')
    title('Right singular vector')
    grid on
end

if run_grid_benchmark
    benchmark_table = benchmark_svd_methods(params,benchmark_Nx_list,benchmark_Ny_list,nsv,opts);
    disp(benchmark_table)
    writetable(benchmark_table,'benchmark_table_continuous_adjoint.csv');

    figure;
    tiledlayout(1,2)

    nexttile
    loglog(benchmark_table.grid_points,benchmark_table.full_total_cpu_s,'ko-','LineWidth',1.5)
    hold on
    loglog(benchmark_table.grid_points,benchmark_table.matrix_free_total_cpu_s,'rs--','LineWidth',1.5)
    legend('full matrix','matrix-free','Location','best')
    xlabel('Nx Ny')
    ylabel('CPU time (s)')
    title('SVD CPU timing')
    grid on

    nexttile
    loglog(benchmark_table.grid_points,benchmark_table.full_memory_MB,'ko-','LineWidth',1.5)
    hold on
    loglog(benchmark_table.grid_points,benchmark_table.matrix_free_memory_MB,'rs--','LineWidth',1.5)
    legend('full matrix','matrix-free','Location','best')
    xlabel('Nx Ny')
    ylabel('memory estimate (MB)')
    title('Dominant array memory')
    grid on
end

if full_matrix_only

    %% ================= FULL MATRIX ONLY =================
    full_matrix_total_cpu_start = cputime;

    % Build full resolvent matrix
    full_matrix_build_cpu_start = cputime;
    [H,L] = resolvent_2D_xy(params);
    full_matrix_build_cpu_time_s = ...
        cputime - full_matrix_build_cpu_start;

    % Compute leading singular value/vector
    full_matrix_svd_cpu_start = cputime;
    [U_full,S_full,V_full] = svds(H,nsv,'largest');
    full_matrix_svd_cpu_time_s = ...
        cputime - full_matrix_svd_cpu_start;

    full_matrix_total_cpu_time_s = ...
        cputime - full_matrix_total_cpu_start;

    sigma_full = diag(S_full);

    %% Display singular value
    disp(table(sigma_full, ...
        'VariableNames',{'sigma_full'}))

    %% CPU timing table
    cpu_time_table = table( ...
        full_matrix_build_cpu_time_s, ...
        full_matrix_svd_cpu_time_s, ...
        full_matrix_total_cpu_time_s, ...
        'VariableNames',{ ...
        'full_matrix_build_cpu_s', ...
        'full_matrix_svds_cpu_s', ...
        'full_matrix_total_cpu_s'});

    disp(cpu_time_table)

end


function benchmark_table = benchmark_svd_methods(params_template,Nx_list,Ny_list,nsv,svds_opts)
if numel(Nx_list) ~= numel(Ny_list)
    error('benchmark_svd_methods:gridListSizeMismatch', ...
        'Nx_list and Ny_list must have the same number of entries.');
end

Nx_cases = Nx_list(:);
Ny_cases = Ny_list(:);
nCases = numel(Nx_cases);

grid_points = Nx_cases.*Ny_cases;
state_size = 3*grid_points;
linear_system_size = 4*grid_points;
full_build_cpu_s = zeros(nCases,1);
full_svds_cpu_s = zeros(nCases,1);
full_total_cpu_s = zeros(nCases,1);
matrix_free_setup_cpu_s = zeros(nCases,1);
matrix_free_svds_cpu_s = zeros(nCases,1);
matrix_free_total_cpu_s = zeros(nCases,1);
sigma_full = zeros(nCases,1);
sigma_matrix_free = zeros(nCases,1);
left_vector_error = zeros(nCases,1);
right_vector_error = zeros(nCases,1);
full_memory_MB = zeros(nCases,1);
matrix_free_memory_MB = zeros(nCases,1);

for iCase = 1:nCases
%     params_i = params_template;
%     params_i.Nx = Nx_cases(iCase);
%     params_i.Ny = Ny_cases(iCase);
%     params_i.N = params_i.Nx*params_i.Ny;
%     params_i = prepare_params(params_i);
    params_i = rebuild_benchmark_params( ...
    params_template,Nx_cases(iCase),Ny_cases(iCase));
    op_size = [3*params_i.Nx*params_i.Ny,3*params_i.Nx*params_i.Ny];

    drawnow;
    full_total_cpu_start = cputime;
    full_build_cpu_start = cputime;
    [H_i,L_i] = resolvent_2D_xy(params_i);
    full_build_cpu_s(iCase) = cputime-full_build_cpu_start;

    full_svds_cpu_start = cputime;
    [U_full_i,S_full_i,V_full_i] = svds(H_i,nsv,'largest',svds_opts);
    full_svds_cpu_s(iCase) = cputime-full_svds_cpu_start;
    full_total_cpu_s(iCase) = cputime-full_total_cpu_start;
    sigma_full(iCase) = S_full_i(1,1);

    full_memory_MB(iCase) = workspace_bytes( ...
        H_i,L_i,U_full_i,S_full_i,V_full_i)/1024^2;

    matrix_free_total_cpu_start = cputime;
    matrix_free_setup_cpu_start = cputime;
    params_i = add_laplacian_preconditioner(params_i);
    matrix_free_setup_cpu_s(iCase) = cputime-matrix_free_setup_cpu_start;

    %clear H_i L_i
    drawnow;

    H_mf_i = @(f,tflag) H_fun(f,tflag,params_i);
    matrix_free_svds_cpu_start = cputime;
    [U_mf_i,S_mf_i,V_mf_i] = svds(H_mf_i,op_size,nsv,'largest',svds_opts);
    matrix_free_svds_cpu_s(iCase) = cputime-matrix_free_svds_cpu_start;
    matrix_free_total_cpu_s(iCase) = cputime-matrix_free_total_cpu_start;
    sigma_matrix_free(iCase) = S_mf_i(1,1);

    [~,u_err_i] = align_singular_vectors(U_mf_i,U_full_i);
    [~,v_err_i] = align_singular_vectors(V_mf_i,V_full_i);

    % nsv = 1, so store the leading-vector errors.
    left_vector_error(iCase) = u_err_i(1);
    right_vector_error(iCase) = v_err_i(1);

    matrix_free_memory_MB(iCase) = (workspace_bytes( ...
        params_i,U_mf_i,S_mf_i,V_mf_i)+matrix_free_workspace_bytes(params_i))/1024^2;

    if Nx_cases(iCase) == Nx_list(end) && Ny_cases(iCase) == Ny_list(end)
        save(sprintf('resolvent_results_mf_%d.mat', Nx_cases(iCase)), ...
            'U_mf_i', 'S_mf_i', 'V_mf_i', 'sigma_matrix_free', 'params_i', ...
            'matrix_free_setup_cpu_s', 'matrix_free_svds_cpu_s', 'matrix_free_total_cpu_s', '-v7.3');
        save(sprintf('resolvent_results_full_%d.mat', Nx_cases(iCase)), ...
            'U_full_i', 'S_full_i', 'V_full_i', 'sigma_full', 'params_i','-v7.3');
    end

    %clear H_mf_i U_mf_i S_mf_i V_mf_i L_i
    clear H_mf_i U_mf_i S_mf_i V_mf_i ...
        U_full_i V_full_i u_err_i v_err_i

end

sigma_relative_error = abs(sigma_matrix_free-sigma_full)./abs(sigma_full);

benchmark_table = table(Nx_cases,Ny_cases,grid_points,state_size, ...
    linear_system_size,full_build_cpu_s,full_svds_cpu_s,full_total_cpu_s, ...
    matrix_free_setup_cpu_s,matrix_free_svds_cpu_s,matrix_free_total_cpu_s, ...
    full_memory_MB,matrix_free_memory_MB,sigma_full,sigma_matrix_free, ...
    sigma_relative_error,left_vector_error,right_vector_error, ...
    'VariableNames',{'Nx','Ny','grid_points','state_size', ...
    'linear_system_size','full_build_cpu_s','full_svds_cpu_s', ...
    'full_total_cpu_s','matrix_free_setup_cpu_s', ...
    'matrix_free_svds_cpu_s','matrix_free_total_cpu_s', ...
    'full_memory_MB','matrix_free_memory_MB','sigma_full', ...
    'sigma_matrix_free','sigma_relative_error', ...
    'left_vector_error','right_vector_error'});
end


function bytes = matrix_free_workspace_bytes(params)
n = 4*params.N; %the size of state vector.
restart = params.gmres_restart;

% Complex double Krylov basis used by GMRES, plus a few work vectors.
bytes_per_complex_double = 16;
bytes = bytes_per_complex_double*(n*(restart+1)+6*n);
if strcmp(params.preconditioner_type,'fft_1d_laplacian')
    bytes = bytes+params.Nx*params.Ny^2*(16+8);
elseif strcmp(params.preconditioner_type,'fft_1d_resolvent')
    % Forward and continuous-adjoint complex factor arrays.
    bytes = bytes+2*params.Nx*(4*params.Ny)^2*16;
end
end


function params=prepare_params(params)
Ny=params.Ny;
Nx=params.Nx;
%     [params.cheb_y,DM]=chebdif(Ny,2);
%     params.Dy_mat=(2/params.Ly)*DM(:,:,1);
%     params.Dyy_mat=(2/params.Ly)^2*DM(:,:,2);
%
%     [x,Dx_mat]=fourdif(Nx,1);
%     [~,Dxx_mat]=fourdif(Nx,2);
%     params.Dx_mat=(2*pi/params.Lx)*Dx_mat;
%     params.Dxx_mat=(2*pi/params.Lx)^2*Dxx_mat;
N=params.N;

%vectorized velocity and velocity gradient that will be used to form
%resolvent operators

% Matrix views used by L_fun.  Keep the vector views above for the
% explicitly assembled full-matrix operator.
% params.U_mat=U_xy_coarse;
% params.V_mat=V_xy_coarse;
% params.dUdx_mat=dUdx_xy_coarse;
% params.dUdy_mat=dUdy_xy_coarse;
% params.dVdx_mat=dVdx_xy_coarse;
% params.dVdy_mat=dVdy_xy_coarse;
% params.K_inv_mat=mask_smooth;
% 
% params.U_xave_diag=U_xave_diag;
% params.V_xave_diag=V_xave_diag;
% params.dUdx_xave_diag=dUdx_xave_diag;
% params.dUdy_xave_diag=dUdy_xave_diag;
% params.dVdx_xave_diag=dVdx_xave_diag;
% params.dVdy_xave_diag=dVdy_xave_diag;
% params.K_inv_xave_diag=K_inv_xave_diag;

[~,w]=clencurt(Ny-1);
params.w=w(:);
params.w_2D=reshape(params.w*ones(1,Nx),N,1);
%% This is the C block matrix . as you can see no pressure is in the 4th row
params.w_all=[params.w_2D;params.w_2D;params.w_2D]; 

end

function params = add_laplacian_preconditioner(params)
params.preconditioner_type = validatestring(params.preconditioner_type, ...
    {'fft_1d_laplacian','fft_1d_resolvent'});

Ny=params.Ny;
Nx=params.Nx;
wave=fourier_wavenumbers(Nx);

if strcmp(params.preconditioner_type,'fft_1d_laplacian')
    Dyy=params.Dyy_mat;
    params.laplacian_1d_factors=cell(Nx,1);
    for kx_ind=1:Nx
        dxx_eigenvalue=-(2*pi/params.Lx)^2*wave(kx_ind)^2;
        lap_k=Dyy+(dxx_eigenvalue-params.kz^2)*speye(Ny);
        lap_k(1,:)=[1,zeros(1,Ny-1)];
        % lap_k(1,1)=1;
        lap_k(Ny,:)=[zeros(1,Ny-1),1];
        % lap_k(Ny,Ny)=1;

        %preconditioning by directly inverse.
        params.laplacian_1d_factors{kx_ind}=inv(lap_k);

        %preconditioning by LU decomposition
        %params.laplacian_1d_factors{kx_ind}=decomposition(lap_k,'lu');
    end
elseif strcmp(params.preconditioner_type,'fft_1d_resolvent')
    block_size=4*Ny;
    params.resolvent_1d_factors=complex(zeros(block_size,block_size,Nx));
    params.resolvent_1d_adjoint_factors=complex(zeros(block_size,block_size,Nx));
    for kx_ind=1:Nx
        resolvent_k=build_1d_resolvent_preconditioner_block(params,wave(kx_ind));
        adjoint_resolvent_k=build_1d_adjoint_resolvent_preconditioner_block(params,wave(kx_ind));

        params.resolvent_1d_factors(:,:,kx_ind)=inv(resolvent_k);
        params.resolvent_1d_adjoint_factors(:,:,kx_ind)=inv(adjoint_resolvent_k);
    end
end
end

function P=build_1d_resolvent_preconditioner_block(params,wave_k)
Ny=params.Ny;
kz=params.kz;
Re=params.Re;
I=speye(Ny);
Z=sparse(Ny,Ny);

kx=(2*pi/params.Lx)*wave_k;
if mod(params.Nx,2)==0 && wave_k==-params.Nx/2
    Dx_k=Z;
else
    Dx_k=1i*kx*I;
end
lap_k=params.Dyy_mat-(kx^2+kz^2)*I;

% This FFT preconditioner is exact in x only for x-homogeneous base
% flow. For x-varying coefficients, use their streamwise average.

U=params.U_xave_diag;
V=params.V_xave_diag;
dUdx=params.dUdx_xave_diag;
dUdy=params.dUdy_xave_diag;
dVdx=params.dVdx_xave_diag;
dVdy=params.dVdy_xave_diag;
K_inv=params.K_inv_xave_diag;

%  fprintf('dUdx diagonal entries:\n');
%     disp(max(mean(dUdx,1)));
% 
%     fprintf('dUdy diagonal entries:\n');
%     disp(max(mean(dUdy,1)));
%     fprintf('dVdx diagonal entries:\n');
%     disp(max(mean(dVdx,1)));
% 
%     fprintf('dVdy diagonal entries:\n');
%     disp(max(mean(dVdy,1)));
%     keyboard

scalar_A=1i*params.omega*I+U*Dx_k+V*params.Dy_mat-lap_k/Re+K_inv;

% P=[scalar_A+dUdx, dUdy, Z, Dx_k; ...
%     dVdx, scalar_A+dVdy, Z, params.Dy_mat; ...
%     Z, Z, scalar_A, 1i*kz*I; ...
%     -Dx_k, -params.Dy_mat, -1i*kz*I, Z];

P=[scalar_A, dUdy, Z, Dx_k; ...
    Z, scalar_A+dVdy, Z, params.Dy_mat; ...
    Z, Z, scalar_A, 1i*kz*I; ...
    -Dx_k, -params.Dy_mat, -1i*kz*I, Z];

wall_rows=[1 Ny];
for block=0:2
    rows=block*Ny+wall_rows;
    P(rows,:)=0;
    P(sub2ind(size(P),rows,rows))=1;
end
end

function P=build_1d_adjoint_resolvent_preconditioner_block(params,wave_k)
Ny=params.Ny;
kz=params.kz;
Re=params.Re;
I=speye(Ny);
Z=sparse(Ny,Ny);

kx=(2*pi/params.Lx)*wave_k;
if mod(params.Nx,2)==0 && wave_k==-params.Nx/2
    Dx_k=Z;
else
    Dx_k=1i*kx*I;
end
lap_k=params.Dyy_mat-(kx^2+kz^2)*I;

% Continuous-adjoint modal operator used by L_fun(...,'transp',...).
% For an x-varying base flow, retain the same streamwise averaging used
% by the forward Fourier-block preconditioner.
U=params.U_xave_diag;
V=params.V_xave_diag;
dUdx=params.dUdx_xave_diag; % changed to zero in this update as streamwise average of grad_x should be 0
dUdy=params.dUdy_xave_diag;
dVdx=params.dVdx_xave_diag; % % changed to zero in this update as streamwise average of grad_x should be 0
dVdy=params.dVdy_xave_diag;
K_inv=params.K_inv_xave_diag;

scalar_A=-1i*params.omega*I-U*Dx_k-V*params.Dy_mat-lap_k/Re+K_inv;

% P=[scalar_A+dUdx, dVdx, Z, Dx_k; ...
%     dUdy, scalar_A+dVdy, Z, params.Dy_mat; ...
%     Z, Z, scalar_A, 1i*kz*I; ...
%     -Dx_k, -params.Dy_mat, -1i*kz*I, Z];

P=[scalar_A, Z, Z, Dx_k; ...
    dUdy, scalar_A+dVdy, Z, params.Dy_mat; ...
    Z, Z, scalar_A, 1i*kz*I; ...
    -Dx_k, -params.Dy_mat, -1i*kz*I, Z];

wall_rows=[1 Ny];
for block=0:2
    rows=block*Ny+wall_rows;
    P(rows,:)=0;
    P(sub2ind(size(P),rows,rows))=1;
end
end


function wave=fourier_wavenumbers(Nx)
N1=floor((Nx-1)/2);
N2=(-Nx/2)*ones(rem(Nx+1,2));
wave=[0:N1 N2 -N1:-1].';
end


function bytes = workspace_bytes(varargin)
bytes = 0;
for k = 1:nargin
    value = varargin{k}; %#ok<NASGU>
    info = whos('value');
    bytes = bytes+info.bytes;
end
end


function [U_aligned,err] = align_singular_vectors(U,U_ref)
U_aligned = U;
err = zeros(1,size(U,2));
for j = 1:size(U,2)
    phase = U_ref(:,j)'*U(:,j);
    if phase ~= 0
        U_aligned(:,j) = U(:,j)*conj(phase)/abs(phase);
    end
    err(j) = norm(U_aligned(:,j)-U_ref(:,j))/norm(U_ref(:,j));
end
end

%% boundary condition and matrix free wrapper for H_mf: returns Hf for forward ('notransp')
%% and H* f for adjoint for 'transp'
%H_f = Hf or H_fun = H*f is the function
function H_f=H_fun(f,tflag,params)

global gmres_info
Ny=params.Ny;
%Nx=params.Nx;
N=params.N;

w_all=params.w_all;

%if strcmp(tflag,'notransp')
Bf=[f.*w_all.^(-1/2);
    zeros(N,1)];
% else
%     Bf=[f.*w_all.^(1/2);
%         zeros(N,1)];
% end

%B.C. for u
left_bc=1:Ny:N;
right_bc=Ny:Ny:N;

Bf(left_bc)=0;
Bf(right_bc)=0;

%B.C. for v
Bf(N+left_bc)=0;
Bf(N+right_bc)=0;

%B.C. for w
Bf(2*N+left_bc)=0;
Bf(2*N+right_bc)=0;

tol = params.gmres_tol;
restart = params.gmres_restart;%min(100,4*Ny);
maxit = params.gmres_maxit;

%% H = CL^-1B ; H* = B*(L*)^-1 C*
%% L_inv_u_p = q , where L_fun x q = Bf; q = [u;v;w;p] so this gives us the input in input-output
%% q = (M^-1L)^-1 M^-1Bf
[L_inv_u_p,flag,relres,iter,resvec] = gmres(@(u_p) L_fun(u_p,tflag,params),Bf, ...
    restart,tol,maxit,@(rhs) laplacian_preconditioner_fun(rhs,tflag,params));

gmres_info(end+1,:) = {height(gmres_info)+1,{tflag}, ...
        iter(1),iter(2),numel(resvec)-1,flag,relres};
% if strcmp(tflag,'notransp')
%% H_f = C*q, C = w_all^(1/2)* diag(I,I,I,0)
%% C = w_all^(1/2)* diag(I,I,I,0)
H_f = w_all.^(1/2).*L_inv_u_p(1:3*N,1); % weighted input q, this is fed into svds
   
% else
%     H_f = w_all.^(-1/2).*L_inv_u_p(1:3*N,1);
% end

end


function z=laplacian_preconditioner_fun(rhs,tflag,params)
%Ny=params.Ny;
%Nx=params.Nx;
N=params.N;

if strcmp(params.preconditioner_type,'fft_1d_laplacian')
    z=zeros(4*N,1);
    z(1:N)=laplacian_1d_fft_solve(rhs(1:N),params);
    z(N+1:2*N)=laplacian_1d_fft_solve(rhs(N+1:2*N),params);
    z(2*N+1:3*N)=laplacian_1d_fft_solve(rhs(2*N+1:3*N),params);
    z(3*N+1:4*N)=rhs(3*N+1:4*N);
elseif strcmp(params.preconditioner_type,'fft_1d_resolvent')
    z=resolvent_1d_fft_solve(rhs,tflag,params);
else
    error('Unknown preconditioner type: %s',params.preconditioner_type)
end

end

function z=laplacian_1d_fft_solve(rhs,params)
Ny=params.Ny;
Nx=params.Nx;

rhs_hat=fft(reshape(rhs,Ny,Nx),[],2);
z_hat=zeros(Ny,Nx);
for kx_ind=1:Nx
    %preconditioning solved based on LU decomposition
    %z_hat(:,kx_ind)=params.laplacian_1d_factors{kx_ind}\rhs_hat(:,kx_ind);

    %preconditioning by directly invert the matrix.
    z_hat(:,kx_ind)=params.laplacian_1d_factors{kx_ind}*rhs_hat(:,kx_ind);
end
z=reshape(ifft(z_hat,[],2),Ny*Nx,1);
end

function z=resolvent_1d_fft_solve(rhs,tflag,params)
Ny=params.Ny;
Nx=params.Nx;

rhs_hat=fft(reshape(rhs,Ny,Nx,4),[],2);
z_hat=zeros(Ny,Nx,4);
if strcmp(tflag,'notransp')
    for kx_ind=1:Nx
        rhs_k=reshape(rhs_hat(:,kx_ind,:),4*Ny,1);
        z_k=params.resolvent_1d_factors(:,:,kx_ind)*rhs_k; % multiplying the 1D resolvent terms of the preconditioner to rhs_k
        z_hat(:,kx_ind,:)=reshape(z_k,Ny,1,4);
    end
else
    for kx_ind=1:Nx
        rhs_k=reshape(rhs_hat(:,kx_ind,:),4*Ny,1);
        z_k=params.resolvent_1d_adjoint_factors(:,:,kx_ind)*rhs_k;
        z_hat(:,kx_ind,:)=reshape(z_k,Ny,1,4);
    end
end
z=reshape(ifft(z_hat,[],2),4*Ny*Nx,1);
end

% L_fun or L_u_p is the A operator as a function being fed to the gmres
function L_u_p=L_fun(u_p,tflag,params)
kz=params.kz;
omega=params.omega;
Ny=params.Ny;
Nx=params.Nx;
Re=params.Re;
N=params.N;

% Reshape each state component once.  All differentiation and operator
% assembly below remain in Ny-by-Nx matrix form.
u=reshape(u_p(1:N),Ny,Nx);
v=reshape(u_p(N+1:2*N),Ny,Nx);
w=reshape(u_p(1+2*N:3*N),Ny,Nx);
p=reshape(u_p(1+3*N:4*N),Ny,Nx);

ux=u*params.Dx_mat.';
vx=v*params.Dx_mat.';
wx=w*params.Dx_mat.';
px=p*params.Dx_mat.';

uy=params.Dy_mat*u;
vy=params.Dy_mat*v;
wy=params.Dy_mat*w;
py=params.Dy_mat*p;

lap_u=params.Dyy_mat*u+u*params.Dxx_mat.'-(kz^2).*u;
lap_v=params.Dyy_mat*v+v*params.Dxx_mat.'-(kz^2).*v;
lap_w=params.Dyy_mat*w+w*params.Dxx_mat.'-(kz^2).*w;

U=params.U_mat;
dUdx=params.dUdx_mat;
dUdy=params.dUdy_mat;
V=params.V_mat;
dVdx=params.dVdx_mat;
dVdy=params.dVdy_mat;
K_inv=params.K_inv_mat;
if strcmp(tflag,'notransp')
    L_u=1i*omega*u+U.*ux+V.*uy-lap_u/Re+K_inv.*u;
    L_v=1i*omega*v+U.*vx+V.*vy-lap_v/Re+K_inv.*v;
    L_w=1i*omega*w+U.*wx+V.*wy-lap_w/Re+K_inv.*w;
    out_u=L_u+dUdx.*u+dUdy.*v+px;
    out_v=L_v+dVdx.*u+dVdy.*v+py;
    out_w=L_w+1i*kz*p;
    out_p=-ux-vy-1i*kz*w;

else
    L_u=-1i*omega*u-U.*ux-V.*uy-lap_u/Re+K_inv.*u;
    L_v=-1i*omega*v-U.*vx-V.*vy-lap_v/Re+K_inv.*v;
    L_w=-1i*omega*w-U.*wx-V.*wy-lap_w/Re+K_inv.*w;
    out_u=L_u+dUdx.*u+dVdx.*v+px;
    out_v=L_v+dUdy.*u+dVdy.*v+py;
    out_w=L_w+1i*kz*p;
    out_p=-ux-vy-1i*kz*w;

end

% Apply the velocity wall rows while the fields are still matrices.
out_u([1 Ny],:)=u([1 Ny],:);
out_v([1 Ny],:)=v([1 Ny],:);
out_w([1 Ny],:)=w([1 Ny],:);

% Vectorize each output only once for GMRES.
L_u_p=[out_u(:);out_v(:);out_w(:);out_p(:)];
end


function [H,L]=resolvent_2D_xy(params)

% f_test = randn(3*N,1) + 1i*randn(3*N,1);
% 
% Hf_mf   = H_fun(f_test,'notransp',params);
% Hf_full = H*f_test;
% 
% H_action_error = norm(Hf_mf-Hf_full)/norm(Hf_full);
% 
% fprintf('Forward H-action error = %.3e\n',H_action_error);

% kx=params.kx;
kz=params.kz;
omega=params.omega;
Ny=params.Ny;
Nx=params.Nx;
Ly=params.Ly;
Lx=params.Lx;

N=Nx*Ny;
Re=params.Re;

Ix=speye(Nx);
Iy=speye(Ny);

% [~,DM] = chebdif(Ny,2);
% D1 = (2/Ly)*DM(1:Ny,1:Ny,1);
% D2 = (2/Ly)^2*DM(1:Ny,1:Ny,2);
% 
% %Dy_2D = kron(Ix,Dy);
% %Dyy_2D = kron(Ix,Dyy);
% 
% [~,Dx] = fourdif(Nx,1);
% Dx = (2*pi/Lx)*Dx;
% 
% [~,Dxx] = fourdif(Nx,2);
% Dxx = (2*pi/Lx)^2*Dxx;
% %Dx_2D = kron(Dx,Iy);
% %Dyy_2D = kron(Dxx,Iy);

% Use the identical discrete operators as the matrix-free formulation.
D1  = params.Dy_mat;
D2  = params.Dyy_mat;
Dx  = params.Dx_mat;
Dxx = params.Dxx_mat;
I = speye(N);

lplc = -kz^2*I + kron(Ix,D2) + kron(Dxx,Iy); % sparse

% U0 = 1-y.^2;
% U1 = -2*y;

U = spdiags(params.U,0,N,N);
dUdx = spdiags(params.dUdx,0,N,N);
dUdy = spdiags(params.dUdy,0,N,N);

V = spdiags(params.V,0,N,N);
dVdx = spdiags(params.dVdx,0,N,N);
dVdy = spdiags(params.dVdy,0,N,N);
K_inv = spdiags(params.K_inv(:),0,N,N);

L11 = 1i*omega*I+U*kron(Dx,Iy) + V*kron(Ix,D1) + dUdx - (1/Re).*lplc + K_inv;
L12 = dUdy;
L21 = dVdx;
L14 = kron(Dx,Iy);
L22 = 1i*omega*I+U*kron(Dx,Iy) + V*kron(Ix,D1) + dVdy - (1/Re)*lplc + K_inv;
L24 = kron(Ix,D1);
L33 = 1i*omega*I+U*kron(Dx,Iy) + V*kron(Ix,D1) - (1/Re)*lplc + K_inv;
L34 = 1i*kz*I;
L41 = -kron(Dx,Iy);
L42 = -kron(Ix,D1);
L43 = -1i*kz*kron(Ix,Iy);

% Assemble A (sparse)
L = [L11 L12 sparse(N,N) L14; ...
    L21 L22 sparse(N,N) L24; ...
    sparse(N,N) sparse(N,N) L33 L34; ...
    L41 L42 L43 sparse(N,N)];

% B = spalloc(4*N,3*N, 3*N);
% C = spalloc(3*N,4*N, 3*N);

[~,w] = clencurt(Ny-1);

w_sqrt = sqrt(w);
w_sqrt = spdiags(w_sqrt',0,Ny,Ny);
w_sqrt = kron(Ix, w_sqrt);

w_sqrt_inv = 1./(sqrt(w)); % elementwise sqrt
w_sqrt_inv = spdiags(w_sqrt_inv',0,Ny,Ny);
w_sqrt_inv = kron(Ix, w_sqrt_inv);

zero_matrix = spalloc(N,N, N);
B = [w_sqrt_inv zero_matrix zero_matrix;
    zero_matrix w_sqrt_inv zero_matrix;
    zero_matrix zero_matrix w_sqrt_inv;
    zero_matrix zero_matrix zero_matrix];

C = [w_sqrt zero_matrix zero_matrix zero_matrix;
    zero_matrix w_sqrt zero_matrix zero_matrix;
    zero_matrix zero_matrix w_sqrt zero_matrix];

left_bc=1:Ny:N;
right_bc=Ny:Ny:N;

%B.C. of u
L(left_bc,:)=0;
L(sub2ind(size(L),left_bc,left_bc))=1;

L(right_bc,:)=0;
L(sub2ind(size(L),right_bc,right_bc))=1;

B(left_bc,:)=0;
B(right_bc,:)=0;

%B.C. of v
L(N+left_bc,:)=0;
L(sub2ind(size(L),N+left_bc,N+left_bc))=1;

L(N+right_bc,:)=0;
L(sub2ind(size(L),N+right_bc,N+right_bc))=1;

B(N+left_bc,:)=0;
B(N+right_bc,:)=0;

%B.C. of w
L(2*N+left_bc,:)=0;
L(sub2ind(size(L),2*N+left_bc,2*N+left_bc))=1;

L(2*N+right_bc,:)=0;
L(sub2ind(size(L),2*N+right_bc,2*N+right_bc))=1;

B(2*N+left_bc,:)=0;
B(2*N+right_bc,:)=0;

H = C*(L\B);

end

function params_i = rebuild_benchmark_params(params_template,Nx,Ny)

params_i = params_template;

params_i.Nx = Nx;
params_i.Ny = Ny;
params_i.N  = Nx*Ny;

Lx = params_i.Lx;
Ly = params_i.Ly;

%% Rebuild Fourier operators
[x_fourier,Dxx] = fourdif(Nx,2);
[~,Dx] = fourdif(Nx,1);

params_i.Dx_mat  = Dx*(2*pi/Lx);
params_i.Dxx_mat = Dxx*(2*pi/Lx)^2;

x_new = x_fourier*(Lx/(2*pi));

%% Rebuild Chebyshev operators
[cheb_y,DM] = chebdif(Ny,2);

params_i.Dy_mat  = -DM(:,:,1)*(2/Ly);
params_i.Dyy_mat =  DM(:,:,2)*(2/Ly)^2;

y_new = flip(cheb_y)*(Ly/2);

params_i.x = x_new;
params_i.cheb_y = y_new;

%% Interpolate fine-grid mean flow
source = params_template.benchmark_source;

[X_old,Y_old] = meshgrid(source.x_old,source.y_old);
[X_new,Y_new] = meshgrid(x_new,y_new);

U_new = interp2(X_old,Y_old,source.U_f, ...
    X_new,Y_new,'linear');

V_new = interp2(X_old,Y_old,source.V_f, ...
    X_new,Y_new,'linear');

dUdx_new = interp2(X_old,Y_old,source.Ux_f, ...
    X_new,Y_new,'linear');

dUdy_new = interp2(X_old,Y_old,source.Uy_f, ...
    X_new,Y_new,'linear');

dVdx_new = interp2(X_old,Y_old,source.Vx_f, ...
    X_new,Y_new,'linear');

dVdy_new = interp2(X_old,Y_old,source.Vy_f, ...
    X_new,Y_new,'linear');

%% Load the mask corresponding to the benchmark grid
mask_filename = sprintf( ...
    'mask_smooth_wavywall_Re%d_%dx%d_epsilon_%.2f.mat', ...
    params_i.Re,Ny,Nx,params_i.epsilon);

if ~isfile(mask_filename)
    error('Benchmark mask file not found:\n%s',mask_filename);
end

mask_data = load(mask_filename,'mask_smooth');
K_new = mask_data.mask_smooth;

if ~isequal(size(K_new),[Ny,Nx])
    error(['Mask in %s has size %d x %d, but benchmark ' ...
        'requires %d x %d.'], ...
        mask_filename,size(K_new,1),size(K_new,2),Ny,Nx);
end

%% Matrix fields
params_i.U_mat     = U_new;
params_i.V_mat     = V_new;
params_i.dUdx_mat  = dUdx_new;
params_i.dUdy_mat  = dUdy_new;
params_i.dVdx_mat  = dVdx_new;
params_i.dVdy_mat  = dVdy_new;
params_i.K_inv_mat = K_new;

%% Vector fields for the assembled operator
params_i.U     = U_new(:);
params_i.V     = V_new(:);
params_i.dUdx  = dUdx_new(:);
params_i.dUdy  = dUdy_new(:);
params_i.dVdx  = dVdx_new(:);
params_i.dVdy  = dVdy_new(:);
params_i.K_inv = K_new(:);

%% Streamwise-averaged fields for 1D FFT preconditioner
params_i.U_xave_diag = spdiags( ...
    mean(U_new,2),0,Ny,Ny);

params_i.V_xave_diag = spdiags( ...
    mean(V_new,2),0,Ny,Ny);

params_i.dUdx_xave_diag = spdiags( ...
    mean(dUdx_new,2),0,Ny,Ny);

params_i.dUdy_xave_diag = spdiags( ...
    mean(dUdy_new,2),0,Ny,Ny);

params_i.dVdx_xave_diag = spdiags( ...
    mean(dVdx_new,2),0,Ny,Ny);

params_i.dVdy_xave_diag = spdiags( ...
    mean(dVdy_new,2),0,Ny,Ny);

params_i.K_inv_xave_diag = spdiags( ...
    mean(K_new,2),0,Ny,Ny);

%% Grid-dependent solver controls and weights
params_i.gmres_maxit = params_i.N;

params_i = prepare_params(params_i);

%% Do not retain the fine-grid source arrays in every benchmark case
params_i = rmfield(params_i,'benchmark_source');

%% Safety checks
assert(isequal(size(params_i.Dy_mat),[Ny,Ny]));
assert(isequal(size(params_i.Dyy_mat),[Ny,Ny]));
assert(isequal(size(params_i.Dx_mat),[Nx,Nx]));
assert(isequal(size(params_i.Dxx_mat),[Nx,Nx]));

assert(isequal(size(params_i.U_mat),[Ny,Nx]));
assert(isequal(size(params_i.V_mat),[Ny,Nx]));
assert(isequal(size(params_i.K_inv_mat),[Ny,Nx]));

assert(all(isfinite(params_i.U_mat),'all'));
assert(all(isfinite(params_i.V_mat),'all'));
assert(all(isfinite(params_i.K_inv_mat),'all'));

end
