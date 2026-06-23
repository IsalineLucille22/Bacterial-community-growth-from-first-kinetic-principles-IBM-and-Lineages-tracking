function intersection_area = compute_rectangle_intersection(Seg, height_cell_1, height_cell_2)
    % A1, B1: Endpoints of the center segment of rectangle 1
    % height_cell_1: Height of rectangle 1
    % A2, B2: Endpoints of the center segment of rectangle 2
    % height_cell_2: Height of rectangle 2

    Seg_1_temp = Seg{1}; %[x(1) x(2); y(1) y(2)]
    Seg_2_temp = Seg{2};
    
    % Define the points A, B for the first segment and C, D for the second segment
    %Put it as argument?
    A1 = Seg_1_temp(:, 1);
    B1 =  Seg_1_temp(:, 2);
    A2 = Seg_2_temp(:, 1);
    B2 = Seg_2_temp(:, 2);
    
    % Compute corners of rectangle 1
    rect1 = get_rectangle_corners(A1, B1, height_cell_1);
    rect2 = get_rectangle_corners(A2, B2, height_cell_2);
    
    % % Create polyshapes for both rectangles
    % poly1 = polyshape(rect1(:,1), rect1(:,2));
    % poly2 = polyshape(rect2(:,1), rect2(:,2));
    % 
    % % Compute intersection polygon
    % intersection = intersect(poly1, poly2);
    % 
    % % Compute intersection area
    % intersection_area = area(intersection);

    intersection_area = rotated_rectangle_intersection(rect1, rect2);
end

function corners = get_rectangle_corners(A, B, height)
    % A, B: Endpoints of the center segment
    % h: Height of the rectangle
    
    % Compute perpendicular vector
    dir_vec = B - A; % Direction vector of the segment
    perp_vec = [-dir_vec(2); dir_vec(1)]; % Perpendicular vector
    perp_vec = (height/2)*perp_vec/norm(perp_vec); % Normalize and scale
    
    % Compute corners
    C1 = A + perp_vec; % Top-left corner
    C2 = A - perp_vec; % Bottom-left corner
    C3 = B + perp_vec; % Top-right corner
    C4 = B - perp_vec; % Bottom-right corner
    
    % Return corners as a matrix
    corners = [C1'; C3'; C4'; C2']; %top-left, top-right, bottom right, bottom left
end

function intersection_area = rotated_rectangle_intersection(rect1, rect2)
    % rect1 and rect2 are 4x2 matrices where each row is a corner (x, y)
    % representing the rectangle in counter-clockwise order.

    % Get the intersection points of the edges
    intersection_points = [];
    for i = 1:4
        % Get edges of rect1
        p1 = rect1(i, :);
        p2 = rect1(mod(i, 4) + 1, :); %Following segment
        for j = 1:4
            % Get edges of rect2
            q1 = rect2(j, :);
            q2 = rect2(mod(j, 4) + 1, :);
            % Compute intersection of edges
            intersect_point = line_segment_intersection(p1, p2, q1, q2);
            if ~isempty(intersect_point)
                intersection_points = [intersection_points; intersect_point];
            end
        end
    end

    % Add corners of rect1 inside rect2
    for i = 1:4
        if Ray_Casting_Algorithm(rect1(i, :), rect2) %If corners are inside
            intersection_points = [intersection_points; rect1(i, :)];
        end
    end

    % Add corners of rect2 inside rect1
    for i = 1:4
        if Ray_Casting_Algorithm(rect2(i, :), rect1)
            intersection_points = [intersection_points; rect2(i, :)];
        end
    end

    % Compute the convex hull of intersection points
    if ~isempty(intersection_points)
        intersection_points = unique(intersection_points, 'rows');
        if size(intersection_points, 1) < 3 || rank(intersection_points - mean(intersection_points, 1)) < 2
            intersection_area = 0; % No intersection or degenerate intersection
        else
            k = convhull(intersection_points(:, 1), intersection_points(:, 2));
            intersection_area = polyarea(intersection_points(k, 1), intersection_points(k, 2));
        end
    else
        intersection_area = 0;
    end
end

function intersection_point = line_segment_intersection(p1, p2, q1, q2)
    % Computes the intersection point of two line segments (p1-p2 and q1-q2)
    % p1, p2: endpoints of the first line segment
    % q1, q2: endpoints of the second line segment
    % Returns the intersection point if it exists, otherwise returns [].

    % Line vectors
    r = p2 - p1; %Segment side 1 rectangle 1
    s = q2 - q1; %Segment side 1 rectangle 1

    % Cross product of r and s
    rxs = cross_product_2D(r, s);

    % Cross products for parameterization
    qp = q1 - p1;
    t = cross_product_2D(qp, s)/rxs;
    u = cross_product_2D(qp, r)/rxs;

    % Check if the intersection is within the bounds of the segments
    if t >= 0 && t <= 1 && u >= 0 && u <= 1
        intersection_point = p1 + t * r; % Intersection point
    elseif abs(rxs) < 1e-10 && cross_product_2D(qp, r) < 1e-10
        %The two segments are collinear
        t_0 = qp*r'/(r*r'); %dot(qp, r)
        t_1 = t_0 + s*r'/(r*r');
        if t_0 > t_1
            temp = t_0;
            t_0 = t_1;
            t_1 = temp;
        end
        t_start = max(0, t_0);
        t_end = min(t_1, 1);
        if t_start <= t_end
            intersection_start = p1 + t_start * r;
            intersection_end = p1 + t_end * r;
            intersection_point = [intersection_start; intersection_end];
        else
            intersection_point = []; % No intersection, lines are parallel
        end
    else
        intersection_point = []; % No intersection within the segments
    end

end

function intersection_point = line_segment_intersection_v2(p1, p2, q1, q2)
    % Computes the intersection point of two line segments (p1-p2 and q1-q2)
    % p1, p2: endpoints of the first line segment
    % q1, q2: endpoints of the second line segment
    % Returns the intersection point if it exists, otherwise returns [].

    % Line vectors
    r = p2 - p1; %Segment side 1 rectangle 1
    s = q2 - q1; %Segment side 1 rectangle 1

    % Cross product of r and s
    rxs = cross_product_2D(r, s);

    % Check if the lines are parallel or coincident
    %If coincide? 
    %Implement the case 1
    if abs(rxs) < 1e-10
        intersection_point = []; % No intersection, lines are parallel
        return;
    end

    % Cross products for parameterization
    qp = q1 - p1;
    t = cross_product_2D(qp, s)/rxs;
    u = cross_product_2D(qp, r)/rxs;

    % Check if the intersection is within the bounds of the segments
    if t >= 0 && t <= 1 && u >= 0 && u <= 1
        intersection_point = p1 + t * r; % Intersection point
    else
        intersection_point = []; % No intersection within the segments
    end
end

function is_inside = Ray_Casting_Algorithm(point, polygon)
    % Checks if a point is inside a polygon
    % point: a 1x2 vector [x, y] representing the point
    % polygon: an Nx2 matrix where each row is a vertex of the polygon [x, y]
    % Returns true if the point is inside the polygon, otherwise false.

    % Extract x and y coordinates
    x = point(1);
    y = point(2);
    poly_x = polygon(:, 1);
    poly_y = polygon(:, 2);

    % Initialize variables
    n = 4;%length(poly_x); % Number of vertices in the polygon %Necessarily 4?
    is_inside = false; % Assume the point is outside

    % Loop through each edge of the polygon
    j = n; % Index of the last vertex 
%     Ray-Casting Algorithm 
    for i = 1:n
        % Check if the point is within the y-bounds of the edge
        if (poly_y(i) > y) ~= (poly_y(j) > y) %Should be inferior to one borne and inferior to the other 
            % Compute the intersection of the edge with the horizontal line at y
            x_intersect = (poly_x(j) - poly_x(i))*(y - poly_y(i))/ ...
                          (poly_y(j) - poly_y(i)) + poly_x(i);

            % Check if the intersection is to the right of the point
            if x < x_intersect
                is_inside = ~is_inside; % Toggle the inside state
            end
        end
        j = i; % Move to the next edge
    end
end


function result = cross_product_2D(v1, v2)
    result = v1(1) * v2(2) - v1(2) * v2(1);
end


