function [mask_smooth, mask_smooth_combined, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface, delta_star, delta_floor, delta_opt] = ...
    normalized_mask_optimizedfn(d_perp, Re, mask_const, c, Ny, Ly, Lx, Nx, Hav, alpha, data_y)

% Defaults for optional args:
if nargin < 10 || isempty(Hav),   Hav   = 1.0; end      % half-height (nondim)

% --- grid floor (your current estimate) ---
j      = 0:(Ny-1);

Y2     = 1 * data_y;
dy_wall_vec = diff(Y2)                % Ny-1 spacings (non-uniform)

dy_wall = min(abs(Y2(2:end) - Y2(1:end-1)))
max_dy_wall  = max(abs(Y2(2:end) - Y2(1:end-1)))
dx_wall = Lx / Nx;
h_wall  = min(dx_wall, dy_wall);
k_cells = 50.0;
delta_floor = k_cells * h_wall;
% --- Burns / "optimal" (Brinkman layer) ---
% eta = c/Re^alpha -> delta_opt = 3.113... * sqrt(c)/Re
%delta_opt = 3.11346786 * (sqrt(c) / Re);
delta_opt = 3.11346786 * sqrt(c/Re^(1+alpha))

% --- choose the thickest (single uniform scalar) ---
delta_star = max([delta_opt, delta_floor]);

solid_region     = (d_perp <= -delta_star);
fluid_region     = (d_perp >=  delta_star);
%fluid_region     = (d_perp > 0);

interface_region = ~(solid_region | fluid_region);
combined_region = (fluid_region | interface_region)

% --- debug prints (matching your Python prints) ---
fprintf('dx_wall = %.10g\n', dx_wall);
fprintf('dy_wall = %.10g\n', dy_wall);
fprintf('h_wall = %.10g\n', h_wall);
fprintf('delta_floor (k_cells*h_wall) = %.10g\n', delta_floor);
fprintf('delta_opt      = %.10g\n', delta_opt);
fprintf('delta_star (max of three)    = %.10g\n', delta_star);

% --- final mask ---
mask_smooth = 0.5 * (1 - erf( sqrt(pi) * (d_perp ./ delta_star) )) * mask_const;

% return region flags
mask_smooth_solid = solid_region;
mask_smooth_fluid = fluid_region;
%mask_smooth_fluid = mask_smooth <1e-3;
mask_smooth_interface = interface_region;
mask_smooth_combined = combined_region;

end
