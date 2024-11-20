% Simulation of 8 points rotating circularly around a central point with Voronoi diagram

% Parameters
central_point = [0, 0]; % Central point coordinates (0, 0)
num_points = 8;         % Number of points

% Define each point's parameters (radius, speed, and initial angle)
points_data = struct('radius', [], 'speed', [], 'angle', []);
radii = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]; % Radii of the orbits
speeds = [0.05, 0.15, 0.065, 0.03, 0.07, 0.04, 0.08, 0.01]; % Angular speeds
angles = linspace(0, 2 * pi, num_points + 1); % Initial angles
angles(end) = []; % Remove duplicate at 360 degrees

% Initialize points_data for each point
for i = 1:num_points
    points_data(i).radius = radii(i);
    points_data(i).speed = speeds(i);
    points_data(i).angle = angles(i);
end

% Colors for each point
colors = lines(num_points + 1); % Include central point

% Set up the figure
figure;
hold on;
axis equal;
axis([-6 6 -6 6]); % Adjust axis limits based on max radius
grid on;
title('Circular Motion with Voronoi Diagram');
xlabel('X');
ylabel('Y');

% Draw central point
central_marker = plot(central_point(1), central_point(2), '*', 'MarkerSize', 8, 'LineWidth', 2);

% Initialize point handles
points = gobjects(num_points, 1);
for i = 1:num_points
    points(i) = plot(nan, nan, '.', 'MarkerSize', 8, 'Color', colors(i, :));
end

% Voronoi diagram handle
voronoi_handle = [];

% Simulation loop
num_frames = 360; % Number of frames in the simulation
for frame = 1:num_frames
    % Collect all positions (central point + rotating points)
    x_positions = central_point(1);
    y_positions = central_point(2);
    
    for i = 1:num_points
        % Compute position of the point
        x = central_point(1) + points_data(i).radius * cos(points_data(i).angle);
        y = central_point(2) + points_data(i).radius * sin(points_data(i).angle);
        
        % Update plot data
        set(points(i), 'XData', x, 'YData', y);
        
        % Update the angle for the next frame
        points_data(i).angle = points_data(i).angle + points_data(i).speed;
        
        % Collect positions
        x_positions = [x_positions, x];
        y_positions = [y_positions, y];
    end
    
    % Draw Voronoi diagram
    if ~isempty(voronoi_handle)
        delete(voronoi_handle); % Clear previous Voronoi lines
    end
    %voronoi_handle = voronoi(x_positions, y_positions, 'k');
    voronoi_handle = voronoiDiagram(x_positions, y_positions);
    
    % Pause to create animation effect
    pause(0.05);
end


