function [h, X_io, X_jo, n] = OverlapVal_v2(Seg, d_0, Over, n, X_o, dist_ind, weights, height_cell_1, height_cell_2, Seg_rect, method_number)
Seg_1_temp = Seg{1}; %[x(1) x(2); y(1) y(2)]
Seg_2_temp = Seg{2};

% Define the points A, B for the first segment and C, D for the second segment
A = Seg_1_temp(:, 1);
B =  Seg_1_temp(:, 2);
pt_1 = Seg_2_temp(:, 1);
pt_2 = Seg_2_temp(:, 2);

% Vectors for each segment
AB = B - A;
pt_1_pt_2 = pt_2 - pt_1;
Apt_1 = A - pt_1;
cos_angle = dot(AB,pt_1_pt_2)/(norm(pt_1_pt_2)*norm(AB));
angle_AB_pt_i = acos(cos_angle);
cos_angle_norm = (cos_angle + 1)/2;

denom = dot(AB, AB)*dot(pt_1_pt_2, pt_1_pt_2) - dot(AB, pt_1_pt_2)*dot(AB, pt_1_pt_2);
length_AB = norm(AB); length_pt_1_pt_2 = norm(pt_1_pt_2);
% ratio_length = length_pt_1_pt_2/length_AB;
ratio_length = min(length_AB / length_pt_1_pt_2, length_pt_1_pt_2 / length_AB);

%If denom is zero, the segments are parallel, handle separately
if denom < 1e-08
    proj = zeros(2, 4); t = zeros(4,1); dist = zeros(4,1);
    n_vec = zeros(2, 4);
    int_pt = [A B pt_1 pt_2];
    % If parallel, find closest point from one endpoint to the other segment
    [proj(:, 1), t(1)] = proj_point_to_segment(A, pt_1, pt_2);
    [proj(:, 2), t(2)] = proj_point_to_segment(B, pt_1, pt_2);
    [proj(:, 3), t(3)] = proj_point_to_segment(pt_1, A, B);
    [proj(:, 4), t(4)] = proj_point_to_segment(pt_2, A, B);
%     dist_proj = norm(proj(:, 4) - proj(:, 3)) + 0*d_0/2;
    n_vec(:, 1) = -(A - proj(:,1)); n_vec(:, 2) = -(B - proj(:,2));
    n_vec(:, 3) = pt_1 - proj(:,3); n_vec(:, 4) = pt_2 - proj(:,4);
    dist(1) = norm(n_vec(:, 1)); % Distance from A to segment pts
    dist(2) = norm(n_vec(:, 2)); % Distance from B to segment pts
    dist(3) = norm(n_vec(:, 3)); % Distance from pt_1 to segment AB
    dist(4) = norm(n_vec(:, 4)); % Distance from pt_2 to segment AB
    [min_dist, ind_min] = min([dist(1), dist(2), dist(3), dist(4)]); % Closest distance
    Over = d_0 - min_dist; 
    n = -n_vec(:, ind_min)/min_dist; %Direction vector
    X_io = proj(:, ind_min)*((2 - ind_min)<0) + int_pt(:, ind_min)*((2 - ind_min)>=0);
    X_jo = proj(:, ind_min)*((2 - ind_min)>=0) + int_pt(:, ind_min)*((2 - ind_min)<0);

    dist_proj = max(0, min(norm(proj(:, 4) - proj(:, 3)), norm(proj(:, 2) - proj(:, 1))));
else
    
    s = (dot(AB, pt_1_pt_2)*dot(pt_1_pt_2, Apt_1) - dot(pt_1_pt_2, pt_1_pt_2)*dot(AB, Apt_1))/denom;
    t = (dot(AB, AB)*dot(pt_1_pt_2, Apt_1) - dot(AB, pt_1_pt_2)*dot(AB, Apt_1))/denom;
    
    % Clamp s and t within [0,1] to ensure points lie on the segments
    s = max(0, min(1, s));
    t = max(0, min(1, t));

    
    % Closest points on each segment
    X_io = A + s*AB;
    X_jo = pt_1 + t*pt_1_pt_2;
    
    min_dist = norm(X_io - X_jo);
    
    % Additional checks for endpoint projections
    % Project endpoints of AB onto CD and vice versa, and calculate distances
    [proj_A_on_CD, ~, dist_A_on_CD] = proj_point_to_segment(A, pt_1, pt_2);
    [proj_B_on_CD, ~, dist_B_on_CD] = proj_point_to_segment(B, pt_1, pt_2);
    [proj_C_on_AB, ~, dist_C_on_AB] = proj_point_to_segment(pt_1, A, B);
    [proj_D_on_AB, ~, dist_D_on_AB] = proj_point_to_segment(pt_2, A, B);
    
    dist_proj = min(norm(proj_D_on_AB - proj_C_on_AB), norm(proj_B_on_CD - proj_A_on_CD)) + 0*d_0/2;

    
    % Find the minimal distance from all cases
    [min_dist, id_case] = min([min_dist, dist_A_on_CD, dist_B_on_CD, dist_C_on_AB, dist_D_on_AB]);
    
    % Update closest points based on the minimal distance case
    switch id_case
    case 2
        X_io = A;
        X_jo = proj_A_on_CD;
    case 3
        X_io = B;
        X_jo = proj_B_on_CD;
    case 4
        X_io = proj_C_on_AB;
        X_jo = pt_1;
    case 5
        X_io = proj_D_on_AB;
        X_jo = pt_2;
    end

    n_vect = X_jo - X_io; %Inverse that
    n = -n_vect/min_dist;
    Over = d_0 - min_dist;
    dist_proj = (dist_proj > 0)*min(max(Over/tan(angle_AB_pt_i), 0), dist_proj);%(dist_proj > 0)*min(max(Over/tan(angle_AB_pt_i), 0), norm(proj_D_on_AB - proj_C_on_AB) + 0*d_0/2);%max(Over/tan(angle_AB_pt_i), 0);%

end
h = max(0, Over);

if method_number == 0
    intersection_area = 0;
    if h > 0
        intersection_area = compute_rectangle_intersection(Seg_rect, height_cell_1, height_cell_2);
    end
    h = intersection_area;
    h = intersection_area + min(intersection_area, h*dist_proj);%h + min(h, h*dist_proj);%h + h*dist_proj + 0*max(0, d_0 - dist_ind);%max(0, h + min(h, Eff_Over)); %max(0, h + h*dist_proj); %max(0, 1/2*h*dist_proj);%max(0, h); %
else
    h = h + min(h, h*dist_proj);
end
if h > 0
    if isnan(n(1))
        n_line = (B - A)/norm(B - A);
        n = [-n_line(2); n_line(1)];
        n = n/norm(n);
    end
else
    n = [1;0];
    X_io = [0;0];
    X_jo = [0;0];
end
end

% Function to project a point onto a segment
function [proj_point, t, min_dist] = proj_point_to_segment(pt, pt_seg_A, pt_seg_B)
    AB = (pt_seg_B - pt_seg_A);
    v_1 = pt - pt_seg_A;
    t = dot(v_1, AB)/dot(AB, AB);
    t = max(0, min(1, t)); % Clamp t to the range [0, 1] to stay on the segment
    proj_point = pt_seg_A + t*AB;
    min_dist = norm(proj_point - pt);
end

% function result = cross_product_2D(v1, v2)
%     result = v1(1) * v2(2) - v1(2) * v2(1);
% end
