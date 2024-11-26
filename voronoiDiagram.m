% Cálculo del diagrama.
% Entrada: posiciones de los puntos (x, y) y el límite del espacio.
% Funcionamiento: revisa cada uno de los puntos y genera el área de Voronoi
% específico de ese punto, almacenandolo en un arreglo de celdas.
% Salida: gráfica del diagrama en el estado actual.
function voronoi_handle = voronoiDiagram(x_positions, y_positions, mass, box)

    S = [x_positions(:), y_positions(:)]; % Une las coordenadas en un arreglo

    cells = cell(size(S, 1), 1);
    
    for p = 1:size(S, 1)
        cell_pts = box;  % El área inicial es el espacio disponible
        Sp = S(p, :); % Punto actual
        for q = 1:size(S, 1)
            if p ~= q
                Sq = S(q, :); % Punto con el que se compara

                [a, b, c] = getPerpendicularBisector(Sp(1), Sp(2), mass(p), Sq(1), Sq(2), mass(q)); % Obtiene los coeficientes del bisector entre los puntos

                [cell_pts, valid] = clipPolygon(cell_pts, a, b, c); % Delimita la nueva área de Voronoi
                if ~valid
                    break;
                end
            end
        end
        cells{p} = [cell_pts; cell_pts(1, :)];
    end

    voronoi_handle = plotVoronoiCells(cells);

    % Calcula los coeficientes de la perpendicular en una posición entre
    % los puntos p y q tal que la fuerza en dicho punto sea 0.
    % Entrada: posiciones y masa de los puntos p y q (x, y, m indica si es
    % coordenada x, y o la masa, respectivamente).
    % Salida: coeficientes de la ecuación del bisector perpendicular.
    function [a, b, c] = getPerpendicularBisector(px, py, pm, qx, qy, qm)
        dis = calculate_distance([px, py], [qx, qy]);

        d1 = dis/(1 + sqrt(qm/pm)); % Distancia respecto al punto p a la que está la zona de fuerza 0.
        
        alfa = calculate_angle([px, py], [qx, qy]);

        x0 = 0;
        y0 = 0;

        % Si p está a la izquierda, se utiliza la distancia directamente y se suma p.
        if(px < qx) 
            x0 = (d1) * cos(alfa) + px;
            y0 = (d1) * sin(alfa) + py;
        % Si q está a la izquierda, el ángulo es referente a este respecto a p, por lo que se toma como punto de partida q,
        % En este caso, la distancia es la complementaria a la de p, es 
        % decir, la distancia total menos la calculada previamente.
        else    
            x0 = (dis-d1) * cos(alfa)+qx;
            y0 = (dis-d1) * sin(alfa)+qy;
        end
        
        %x0 = (px + qx) / 2;
        %y0 = (py + qy) / 2;

        a = qx - px;
        b = qy - py;
        c = -(a * x0 + b * y0);
        
    end

    % De corresponder al área de Voronoi, agrega el bisector dentro del
    % poliígono que representa dicha área.
    % Entrada: polígono actual y los coeficientes del bisector nuevo.
    % Salida: arreglo de puntos que delimitan el polígono y la
    % especificación si es o no válido.
    function [clipped_poly, valid] = clipPolygon(poly, a, b, c)
        distances = a * poly(:, 1) + b * poly(:, 2) + c;
        inside = distances <= 0; % Puntos con distancia negativa están dentro del área de Voronoi
        clipped_poly = [];
        num_vertices = size(poly, 1);
        for i = 1:num_vertices  % Se itera por todos los vértices para determinar si aun son válidos, no lo son, o se deben de adaptar.
            curr = poly(i, :);  % Obtiene el punto actual
            next = poly(mod(i, num_vertices) + 1, :); % Obtiene el punto siguiente (lo considera con módulo para cerrar la figura cuando se llegue al final)
            curr_in = inside(i);    % Indica si el punto aún es válido o no.
            next_in = inside(mod(i, num_vertices) + 1); 

            if curr_in  % Si el punto es válido, lo añade.
                clipped_poly(end + 1, :) = curr; %#ok<AGROW>
            end

            if curr_in ~= next_in % Si entre dos puntos continuos uno es válido y el otro no, se actualiza.
                t = distances(i) / (distances(i) - distances(mod(i, num_vertices) + 1));
                intersection = curr + t * (next - curr);
                clipped_poly(end + 1, :) = intersection; %#ok<AGROW>
            end
        end

        valid = true;
    end
    
    % Genera la gráfica.
    % Entrada: celdas (áreas de Voronoi)
    % Salida: objeto que contiene la gráfica.
    function voronoi_handle = plotVoronoiCells(cells)
        hold on;
        colors = lines(numel(cells));
        voronoi_handle = gobjects(numel(cells), 1);
        for i = 1:numel(cells)
            if ~isempty(cells{i})
                voronoi_handle(i) = plot(cells{i}(:, 1), cells{i}(:, 2), 'Color', colors(i, :), 'LineWidth', 1); % Agrega el área de Voronoi a la gráfica
            end
        end
        hold off;
    end
    
    % Calcula la distancia entre dos puntos.
    % Entrada: arreglos de coordenadas de la forma [x, y] de los puntos a
    % comparar.
    % Salida: distancia en unidades del gráfico.
    function dis = calculate_distance(p1,p2)
        dis = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
    end

    % Calcula el ángulo del punto más a la izquierda respecto al punto más a la derecha.
    % Entrada: arreglos de coordenadas de la forma [x, y] de los puntos a
    % comparar.
    % Salida: ángulo en radianes.
    function [alfa] = calculate_angle(p1,p2)
        alfa = atan((p2(2)-p1(2)) / (p2(1) - p1(1))); 
    end

end
