
central_point = [0, 0];
padding = 0.5;

aux =                   input("Posición x del centro                   : ");
central_point(1) = aux;
aux =                   input("Posición y del centro                   : ");
central_point(1) = aux;

num_points = input("Número de puntos : ");

radii = [];
speeds = [];
angles = [];

for i=1 : num_points
    radii(i) =  input("Indica la distancia del centro           : ");
    speeds(i) = input("Indica la velocidad angular (°/frame)    : ");
end
angles = linspace(0, 2 * pi, num_points + 1);
angles(end) = [];

points_data = struct('radius', num2cell(radii), 'speed', num2cell(speeds), 'angle', num2cell(angles));

% Colores, incluye el punto central
colors = lines(num_points + 1);

min_x = -max(radii);
min_y = -max(radii);
max_x = max(radii);
max_y = max(radii);
box = [min_x - padding, min_y - padding;
           max_x + padding, min_y - padding;
           max_x + padding, max_y + padding;
           min_x - padding, max_y + padding];

figure;
hold on;
axis equal;
axis([ min_x - padding max_x+padding min_y-padding max_y+padding]);
grid on;
title('Diagrama de Voronoi Dinamico');
xlabel('X');
ylabel('Y');

central_marker = plot(central_point(1), central_point(2), '*', 'MarkerSize', 8, 'LineWidth', 2);
points = arrayfun(@(i) plot(nan, nan, '.', 'MarkerSize', 8, 'Color', colors(i, :)), 1:num_points);

voronoi_handle = [];

num_frames = 360;

for frame = 1:num_frames
    x_positions = zeros(num_points + 1, 1);
    y_positions = zeros(num_points + 1, 1);
    x_positions(1) = central_point(1);
    y_positions(1) = central_point(2);

    for i = 1:num_points
        x = central_point(1) + points_data(i).radius * cos(points_data(i).angle);
        y = central_point(2) + points_data(i).radius * sin(points_data(i).angle);

        set(points(i), 'XData', x, 'YData', y);

        points_data(i).angle = points_data(i).angle + points_data(i).speed;

        x_positions(i + 1) = x;
        y_positions(i + 1) = y;
    end

    if ~isempty(voronoi_handle)
        delete(voronoi_handle);
    end
    voronoi_handle = voronoiDiagram(x_positions, y_positions, box);
    drawnow;
end
