function voronoi_handle = voronoiDiagram(x_positions, y_positions)
    % Combine x_positions and y_positions into S
    S = [x_positions(:), y_positions(:)];
    % Define the array S
    % S = [1, 3; 10, 5; 1, 4; 2, 5; 7, 10; 9, 8];
    % S = [1, 3; 1, 4; 2, 5; 7, 10; 9, 8; 10, 5];
    % S = [1, 3; 10, 5; 5, 6];
    % S = [1, 3; 5, 6; 10, 5]
    
    % Define padding
    padding = 0.5;
    
    % Get minimum and maximum values for x and y columns
    min_x = min(S(:, 1));
    min_y = min(S(:, 2));
    max_x = max(S(:, 1));
    max_y = max(S(:, 2));
    
    % Define the box coordinates (clockwise order)
    box = [min_x - padding, min_y - padding;
           max_x + padding, min_y - padding;
           max_x + padding, max_y + padding;
           min_x - padding, max_y + padding];
    
    % Plot the box
    figure;
    hold on;
    
    % Draw the box
    plot([box(:,1); box(1,1)], [box(:,2); box(1,2)], 'k', 'LineWidth', 1.5);
    
    % Plot the original points inside the box
    scatter(S(:, 1), S(:, 2), '.');
    
    % Add grid and labels for better visualization
    grid on;
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Voronoi diagram');
    xlim([min_x - padding, max_x + padding]);
    ylim([min_y - padding, max_y + padding]);
    
    cells = [];
    
    % Colors for plotting
    colors = lines(size(S,1));
    
    % Loop through all points to calculate Voronoi cells
    for p = 1:size(S, 1)
        cell_pts = box;  % Start with the bounding box for each point
        for q = 1:size(S, 1)
            if p ~= q
                intersections = [];  % Initialize intersections array
    
                % Get perpendicular bisector coefficients (a, b, c)
                [a, b, c] = getPerpendicularBisector(S(p, 1), S(p, 2), S(q, 1), S(q, 2));
                N = size(cell_pts, 1);  % Number of vertices in the current polygon
                for i = 1:N
                    % Current edge is from cell_pts(i) to cell_pts(i+1)
                    p1 = cell_pts(i, :);
                    p2 = cell_pts(mod(i, N) + 1, :); % Wrap around
    
                    % Find intersection of the edge with the line
                    [intersection, x, y] = segment_line_intersection(p1, p2, a, b, c);
    
                    if intersection
                        % Append intersection as a row [x, y, edge_index]
                        intersections = [intersections; x, y, i];
                    end
                end
    
                if size(intersections, 1) == 2
    
                    % Take the first two intersections
                    t1 = intersections(1, 1:2);
                    t2 = intersections(2, 1:2);
                    idx1 = intersections(1, 3);
                    idx2 = intersections(2, 3);
    
                    % Construct new_cell
                    N = size(cell_pts, 1);
                    if idx1 < idx2
                        idx_between = idx1+1:idx2;
                    else
                        idx_between = [idx1+1:N, 1:idx2];
                    end
                    new_cell = [t1; cell_pts(idx_between, :); t2];
    
                    % Check if point p is inside new_cell
                    px = S(p, 1);
                    py = S(p, 2);
                    in = inpolygon(px, py, new_cell(:, 1), new_cell(:, 2));
    
                    if in
                        cell_pts = new_cell;
                    else
                        % Reverse the order
                        if idx2 < idx1
                            idx_between = idx2+1:idx1;
                        else
                            idx_between = [idx2+1:N, 1:idx1];
                        end
                        new_cell = [t2; cell_pts(idx_between, :); t1];
                        cell_pts = new_cell;
                        %cell_pts = [cell_pts; cell_pts(1, :)];
                    end
                end
            end
        end
    
        % Close the polygon
        cell_pts = [cell_pts; cell_pts(1, :)];
    
        % Store the resulting cell
        cells{p} = cell_pts;
    
        % Plot the cell
        plot(cell_pts(:, 1), cell_pts(:, 2), 'Color', colors(p, :), 'LineWidth', 1.5);
    end
    
    hold off;
    voronoi_handle = gcf;
    
    % Function to get perpendicular bisector coefficients
    function [a, b, c] = getPerpendicularBisector(px, py, qx, qy)
        % Midpoint of line segment
        x0 = (px + qx) / 2;
        y0 = (py + qy) / 2;
    
        a = qx - px;
        b = qy - py;
        
        % Line equation: ax + by + c = 0 passing through midpoint
        c = -(a * x0 + b * y0);
    end
    
    % Function to find intersection between segment and line
    function [intersection, x, y] = segment_line_intersection(A, B, a, b, c)
        % Line segment AB
        x1 = A(1); y1 = A(2);
        x2 = B(1); y2 = B(2);
    
        % Line coefficients for segment AB
        denom = a*(x2 - x1) + b*(y2 - y1);
    
        if denom == 0
            % Line and segment are parallel
            intersection = false;
            x = NaN;
            y = NaN;
            return;
        end
    
        t = -(a*x1 + b*y1 + c) / denom;
    
        if t >= 0 && t <= 1
            x = x1 + t*(x2 - x1);
            y = y1 + t*(y2 - y1);
            intersection = true;
        else
            intersection = false;
            x = NaN;
            y = NaN;
        end
    end
end

