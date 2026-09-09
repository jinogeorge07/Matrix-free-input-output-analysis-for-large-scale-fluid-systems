function [d_perp] = compute_signed_distance_mask(x_grid, y_grid, Lx, y0, A1, A2)
    % Input:
    % x_grid, y_grid: 2D meshgrid of x and y coordinates (size: ny x nx)
    % Lx: Domain length in x
    % y0: Wall center (e.g., 0.6)
    % A1, A2: Amplitudes for lower and upper wall
    % mask_threshold: transition width (e.g., 3*dy)

    nsample = 6000;
    s_vals = linspace(0, Lx, nsample);
    y1_vals = -y0 - A1 * sin(2 * pi * s_vals / Lx); % bottom wall
    y2_vals =  y0 + A2 * sin(2 * pi * s_vals / Lx); % top wall

    [ny, nx] = size(x_grid);
    d_perp = zeros(ny, nx);

    for j = 1:ny
        for i = 1:nx
            x_pt = x_grid(j, i);
            y_pt = y_grid(j, i);

            % Distances to bottom wall
            dist2_1 = (x_pt - s_vals).^2 + (y_pt - y1_vals).^2;
            [min_dist1, idx1] = min(dist2_1);
            y1_closest = y1_vals(idx1);
            d1 = sqrt(min_dist1);

            % Distances to top wall
            dist2_2 = (x_pt - s_vals).^2 + (y_pt - y2_vals).^2;
            [min_dist2, idx2] = min(dist2_2);
            y2_closest = y2_vals(idx2);
            d2 = sqrt(min_dist2);

            % Determine signed distance
            if d1 < d2
                sign_val = 1 * (y_pt > y1_closest) - 1 * (y_pt <= y1_closest);
                d_perp(j, i) = sign_val * d1;
            else
                sign_val = 1 * (y_pt < y2_closest) - 1 * (y_pt >= y2_closest);
                d_perp(j, i) = sign_val * d2;
            end
        end
    end

end

