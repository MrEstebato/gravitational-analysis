% Define the array S
S = [1, 3; 10, 5; 1, 4; 2, 5; 7, 10; 9, 8];

% voronoi(S(:, 1), S(:, 2))

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
scatter(S(:, 1), S(:, 2), 'filled');

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
                    % Append intersection as a row [x, y]
                    intersections = [intersections; x, y];

                    % Debug: Plot the intersection point
                    % plot(x, y, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');  % Blue circle for intersection
                    
                    % Break if we have exactly two intersections
                    if size(intersections, 1) == 2
                        break;
                    end
                end
            end
            
            % If we have two intersections, construct a new polygon
            if size(intersections, 1) == 2
                % Extract vertices
                t1 = intersections(1, :);
                t2 = intersections(2, :);

                % Find vertices between t1 and t2 in the current polygon
                idx1 = findPointIndex(cell_pts, t1);
                idx2 = findPointIndex(cell_pts, t2);

                % Check orientation of t1 and t2 in the polygon
                if idx1 < idx2
                    new_cell = [t1; cell_pts(idx1 + 1:idx2, :); t2];
                else
                    new_cell = [t1; cell_pts(idx1 + 1:end, :); cell_pts(1:idx2, :); t2];
                end

                % Check if point p is inside new_cell
                px = S(p, 1);
                py = S(p, 2);
                in = inpolygon(px, py, new_cell(:, 1), new_cell(:, 2));

                if in
                    cell_pts = new_cell;
                else
                    % Reverse the polygon if p is not inside
                    new_cell = [t2; cell_pts(idx2 + 1:end, :); cell_pts(1:idx1, :); t1];
                    cell_pts = new_cell;
                end
            end
        end
    end
    
    % Store the resulting cell
    cells{p} = cell_pts;
    
    % Close the polygon for visualization
    %cell_pts = [cell_pts; cell_pts(1, :)];
    
    % Plot the cell
    fill(cell_pts(:, 1), cell_pts(:, 2), colors(p, :), 'FaceAlpha', 0.2, 'EdgeColor', colors(p, :));
end

hold off;

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

function [intersection, x, y] = segment_line_intersection(A, B, a, b, c)
    % Input: 
    % A = [Ax, Ay] -> First point of segment AB
    % B = [Bx, By] -> Second point of segment AB
    % a, b, c -> Coefficients of the line ax + by + c = 0
    
    % Output:
    % intersection -> Boolean value, true if intersection exists and lies on the segment
    % x, y -> Coordinates of the intersection point (NaN if no intersection)
    
    % Unpack inputs
    Ax = A(1); Ay = A(2);
    Bx = B(1); By = B(2);
    
    % Calculate vector u and v
    u = [a, b, c];
    v = [Ay - By, Bx - Ax, Ax*By - Bx*Ay];
    
    % Compute cross product u x v
    z = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)];
    z1 = z(1); z2 = z(2); z3 = z(3);
    
    % Check if intersection exists
    if z3 == 0
        % Parallel or coincident
        intersection = false;
        x = NaN;
        y = NaN;
        return;
    end
    
    % Intersection point
    x = z1 / z3;
    y = z2 / z3;
    
    % Define bounds for the segment
    x1 = min(Ax, Bx);
    x2 = max(Ax, Bx);
    y1 = min(Ay, By);
    y2 = max(Ay, By);
    
    % Check if the intersection point lies on the segment AB
    if (x1 <= x && x <= x2) && (y1 <= y && y <= y2)
        intersection = true;
    else
        intersection = false;
        x = NaN;
        y = NaN;
    end
end

function idx = findPointIndex(cell_pts, point)
    % Find the index of the vertex in cell_pts closest to 'point'
    distances = sqrt(sum((cell_pts - point).^2, 2));
    [~, idx] = min(distances);
end


