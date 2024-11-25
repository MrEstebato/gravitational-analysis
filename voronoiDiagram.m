function voronoi_handle = voronoiDiagram(x_positions, y_positions, box)

    S = [x_positions(:), y_positions(:)];

    cells = cell(size(S, 1), 1);
    
    for p = 1:size(S, 1)
        cell_pts = box; 
        Sp = S(p, :);
        for q = 1:size(S, 1)
            if p ~= q
                Sq = S(q, :);

                [a, b, c] = getPerpendicularBisector(Sp(1), Sp(2), Sq(1), Sq(2));

                [cell_pts, valid] = clipPolygon(cell_pts, a, b, c);
                if ~valid
                    break;
                end
            end
        end
        cells{p} = [cell_pts; cell_pts(1, :)];
    end

    voronoi_handle = plotVoronoiCells(cells);

    function [a, b, c] = getPerpendicularBisector(px, py, qx, qy)

        x0 = (px + qx) / 2;
        y0 = (py + qy) / 2;

        a = qx - px;
        b = qy - py;
        c = -(a * x0 + b * y0);
        
    end

    function [clipped_poly, valid] = clipPolygon(poly, a, b, c)
        distances = a * poly(:, 1) + b * poly(:, 2) + c;
        inside = distances <= 0;
        clipped_poly = [];
        num_vertices = size(poly, 1);
        for i = 1:num_vertices
            curr = poly(i, :);
            next = poly(mod(i, num_vertices) + 1, :);
            curr_in = inside(i);
            next_in = inside(mod(i, num_vertices) + 1);

            if curr_in
                clipped_poly(end + 1, :) = curr; %#ok<AGROW>
            end

            if curr_in ~= next_in
                t = distances(i) / (distances(i) - distances(mod(i, num_vertices) + 1));
                intersection = curr + t * (next - curr);
                clipped_poly(end + 1, :) = intersection; %#ok<AGROW>
            end
        end
        valid = true;
    end

    function voronoi_handle = plotVoronoiCells(cells)
        hold on;
        colors = lines(numel(cells));
        voronoi_handle = gobjects(numel(cells), 1);
        for i = 1:numel(cells)
            if ~isempty(cells{i})
                voronoi_handle(i) = plot(cells{i}(:, 1), cells{i}(:, 2), 'Color', colors(i, :), 'LineWidth', 1);
            end
        end
        hold off;
    end
end
