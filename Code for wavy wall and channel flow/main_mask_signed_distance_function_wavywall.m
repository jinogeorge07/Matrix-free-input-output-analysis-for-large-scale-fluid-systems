%% main signed distance function , binary mask and normalized mask for input putput SVDS
% Define domain
%Lx = 0.6 * pi;

% %% Read dedalus data file for x y dataset for grid and velocity profile
% %% ---------------------xxxxxxxxxxxxxxxxxxxxxxxxx--------------------------
%folderpath = 'E:/dedalus/Wavy Wall/dedalus_local_Re75_3/snapshots_channel/mean_v_Re75.00_c1.00' 
%folderpath = 'E:\dedalus\Wavy Wall\flexible signed distance function\optimum_vpm_symmetric_wavywall\epsilon 0.12\dedalus_22423155_Re500\snapshots_channel\mean_v_Re500.00_c1.00';
%folderpath = 'E:\dedalus\Wavy Wall\flexible signed distance function\optimum_vpm_symmetric_wavywall\epsilon 0.12\dedalus_45035570_Re500_NxNy_608x640\snapshots_channel\mean_v_Re500.00_c1.00';

% %% ========================= USER INPUTS ==================================
folderpath = './'
cd(folderpath)


load data_x.mat
load data_y.mat

%Ny = 60; Nx = 60;  % Nz to Nx for streamwise terms

% Wavy validation case Choo parameters
epsilon = 0.12;
h = 1;
y0 = h;

A1 = epsilon;
A2 = epsilon;
Ly = h + 1.1*epsilon;

[x,Dx] = fourdif(Nx,1); % first derivative x
[y_cheb,DM] = chebdif(Ny,2);
y = y_cheb*Ly;

x = x * (Lx / (2*pi)); 

[x_grid, y_grid] = meshgrid(x, y);

% Wall parameters
y0 = 1.0;
A1 = epsilon;
A2 = epsilon;
dy = Ly / Ny;
%mask_threshold = 2 * 0.5 * dy; % 8 to 6 to 3
mask_option = 2;

c = 1.0; 
mask_const = Re/c;  % Example value for penalization stiffness: 400
eta       = c*(1.0/Re)
gamma = sqrt(eta / Re);                      % ε = sqrt(η / Re)
mask_threshold = 3.11346786 * gamma;         % δ* = 3.113... * ε

%% Compute distance and mask
[d_perp] = compute_signed_distance_mask(x_grid, y_grid, Lx, y0, A1, A2);

Hav = y0;
alpha = 1.0;
[mask_smooth, mask_smooth_combined, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface, delta_star, delta_floor, delta_opt] = ...
normalized_mask_optimizedfn(d_perp, Re, mask_const, c, Ny, Ly, Lx, Nx, Hav, alpha,data_y)

% if mask_option == 1
%     %% Normalized smooth mask using erf
%     [mask_smooth,mask_smooth_solid] = normalized_mask_erf(d_perp, mask_threshold, mask_const);
% 
% elseif mask_option == 2
%      %% Normalized optimized smooth mask using peclet number and optimized erf
%     Hav = y0;
%     alpha = 1.0;
%     [mask_smooth, mask_smooth_combined, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface, delta_star, delta_floor, delta_opt] = ...
%     normalized_mask_optimizedfn(d_perp, Re, mask_const, c, Ny, Ly, Lx, Nx, Hav, alpha,data_y)
% end

%% Coordinates of solid region from mask_smooth_solid
[x_solid, y_solid] = deal(x_grid(mask_smooth_solid), y_grid(mask_smooth_solid));
[x_interface, y_interface] = deal(x_grid(mask_smooth_interface), y_grid(mask_smooth_interface));

solid_nodes = [x_solid y_solid];

%% Coordinates of solid region from mask_smooth_solid
[x_fluid, y_fluid] = deal(x_grid(mask_smooth_fluid), y_grid(mask_smooth_fluid));

fluid_nodes = [x_fluid y_fluid];


% Plot2
figure(2);
pcolor(x_grid, y_grid, mask_smooth); shading flat; axis equal tight;
colorbar; title(sprintf('Smooth Normalized Mask (erf) at Re = %.2f', Re));
outdir      = folderpath;  % or wherever you like
plot_fname  = sprintf('smooth_mask_Re%.2f_mask%d.png', Re, mask_const);
full_plotfn = fullfile(outdir, plot_fname);
saveas(gcf, full_plotfn);

% plot_filename = sprintf('smooth normalied mask(erf) for Re%.2f and mask %d.png', Re, mask_const);
% saveas(gcf, plot_filename);

% % Plot3
% figure(3);
% pcolor(x_grid, y_grid, mask_smooth_tanh); shading flat; axis equal tight;
% colorbar; title('Smooth Normalized Mask (tanh)');

% Plot4
figure(4)
plot(x_solid,y_solid,'*');
hold on
% Plot5
plot(x_fluid,y_fluid,'o');
title("solid and fluid cells")
legend('solid cells','fluid cells','interface')
xlim([0, Lx]);
ylim([-Ly, Ly]);
hold on

  %dlmwrite('mask_smooth_erf.csv',mask_smooth); 
% dlmwrite('mask_smooth_tanh.csv',mask_smooth_tanh);
% dlmwrite('solid_nodes.csv',solid_nodes); 

% % Output files
%   save("mask_smooth_hpc_120x96.mat", "mask_smooth");
%  % save("mask_smooth_hpc_80x64.mat");

filename = ['mask_smooth_wavywall_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '_epsilon_' num2str(epsilon) '.mat'];
save(filename, 'mask_smooth');

filename = ['mask_smooth_solid_wavywall_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '.mat'];
save(filename, 'mask_smooth_solid');

  % save("mask_smooth_fluid");
  save("mask_smooth_solid.mat");
  save("x_solid.mat","x_solid");
  save("y_solid.mat","y_solid");


function [mask_smooth, mask_smooth_solid] = normalized_mask_erf(d_perp, mask_threshold, mask_const)
    % Compute a smooth penalization mask based on signed distance
    % - Solid: mask_const
    % - Fluid: 0
    % - Interface: erf transition

    steepness = 1;  % Can increase for sharper interface (e.g., 2), or decrease for smoother

    mask_smooth = zeros(size(d_perp));

    % Region masks
    solid_region     = (d_perp <= -mask_threshold);
    fluid_region     = (d_perp >  mask_threshold);
    interface_region = ~solid_region & ~fluid_region;
        
    mask_smooth = 0.5 * (1 - erf(steepness * sqrt(pi) * d_perp(interface_region) / mask_threshold)) * mask_const;
    % Output solid region mask (binary, same shape)
    mask_smooth_solid = solid_region;
    
end


function mask_smooth = normalized_mask_tanh(d_perp, mask_threshold, mask_const)
    % Compute smooth mask using tanh transition (alternative to erf)
    % - Solid: mask_const
    % - Fluid: 0
    % - Interface: smooth transition via tanh

    mask_smooth = zeros(size(d_perp));

    % Region masks
    solid_region     = (d_perp <= -mask_threshold);
    fluid_region     = (d_perp >  mask_threshold);
    interface_region = ~solid_region & ~fluid_region;

    % Assign values
    mask_smooth(solid_region) = mask_const;
    mask_smooth(fluid_region) = 0;
    mask_smooth(interface_region) = ...
        0.5 * (1 - tanh(d_perp(interface_region) / mask_threshold)) * mask_const;
end

function [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = normalized_mask1_shifted(d_perp, mask_threshold, mask_const)
%NORMALIZED_MASK1_SHIFTED   Shifted erf mask + region indicators (O(ε²) error)
%
%   [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = ...
%       normalized_mask1_shifted(d_perp, mask_threshold, mask_const)
%
%   Inputs:
%     d_perp         Array of signed distances from the wall
%     mask_threshold Half‑width δ of the original erf ramp
%     mask_const     Penalty strength (1/τ)
%
%   Outputs:
%     mask_smooth        Shifted erf mask array
%     mask_smooth_solid  Logical array: true inside the solid region
%     mask_smooth_fluid  Logical array: true inside the fluid region

    % 1) compute displacement length ℓ = δ / √π
    shift = mask_threshold / sqrt(pi);

    % 2) shift distance deeper into the solid
    d_s = d_perp + shift;

    % 3) define regions on the shifted coordinate
    solid_region = (d_s <= -mask_threshold);
    fluid_region = (d_s >=  mask_threshold);
    interface_region  = ~(solid_region | fluid_region);

    % 4) build the shifted erf mask
    % 4) smooth mask (not piecewise anymore)
    mask_smooth = 0.5 * (1 - erf(sqrt(pi) * d_s / mask_threshold)) * mask_const;

    mask_smooth_solid = solid_region;
    mask_smooth_fluid = fluid_region;
    %mask_smooth_interface = interface_region;

end
