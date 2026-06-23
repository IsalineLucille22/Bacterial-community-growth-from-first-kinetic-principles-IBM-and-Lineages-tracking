function [Pos_X, cell_angle, ind_Pred, coll_ind] = NewPos(dt_forces, Pos_X, vect_Cell_Length_tot, Pos_Boundaries, mat_Pred, height_cell, Mass_Cell, Vect_Boundaries, nb_boundaries, t, coll_ind)
vect_Cell_length = vect_Cell_Length_tot(1,:); %Length of a cell
sp_nb = vect_Cell_Length_tot(2,:); %Index for the species
cell_angle = vect_Cell_Length_tot(3,:); %Angle value
g = 0;%If 0 no force due to gravity. The gravity acceleration has to be in um/hours^2
%%%%%%%%Boundary conditions
% Pos_X = [Pos_X Pos_Boundaries'];
% vect_Cell_length = [vect_Cell_length, Vect_Boundaries(1,:)];
% cell_angle = [cell_angle, Vect_Boundaries(2,:)];
% height_cell = [height_cell; Vect_Boundaries(3,:)'];
% sp_nb = [sp_nb, Vect_Boundaries(4,:)];
% Mass_Cell = [Mass_Cell 100000 100000 100000 100000];
%%%%%%%%
n = length(Pos_X(1,:)); %Number of cells at time t
coll_ind_temp = [coll_ind, n - 1]; %When collision with bottom?
h_mat_init = zeros(n, 1);
X_io_mat_init = zeros(n, 2);
ind_Pred = zeros(1,n);
Over = zeros(1,2);
n_init = zeros(2,2);
X_o_init = {n_init, n_init};
fric_param = 0.6;
fric_coeff = 300;
fric_coeff_par = fric_param*fric_coeff;
fric_coeff_per = fric_coeff/fric_param; %Increase friction parameter decrease variation of the angle
Seg_tot = arrayfun(@(x) Rect2Seg([Pos_X(1,x) Pos_X(2,x) vect_Cell_length(x) height_cell(sp_nb(x))], cell_angle(x)),1:n,'UniformOutput',false);
[h_mat, X_io_mat_x, X_io_mat_y, n_ij_mat_x, n_ij_mat_y]  = deal(zeros(n, n));
[h_ij, X_io,X_jo, n_ij] = deal(cell(1,n));
F_ij_init = zeros(2,n);%zeros(2,length(h_ij_temp));
T_ij_init = zeros(1,n);%zeros(1,length(h_ij_temp));
for i = 1:n %(n - nb_boundaries + 4)%If no boundary condition, put n otherwise (n-1) or (n-2)
    X_io_grav = Pos_X(:,i);
    [h_ij{i}, X_io{i}, X_jo{i}, n_ij{i}, X_io_grav] = HertzForce(Pos_X, i, vect_Cell_length, cell_angle, height_cell, Seg_tot, sp_nb, n, h_mat_init, X_io_mat_init, Over, n_init, X_o_init, X_io_grav);
    h_ij_temp = h_ij{i};
    h_mat((i+1):n, i) = h_ij_temp((i+1):n);
    h_mat(i, (i+1):n) = h_mat((i+1):n, i);
    X_io_mat_temp = X_io{i}; X_jo_mat_temp = X_jo{i};
    X_io_mat_x((i+1):n, i) = X_io_mat_temp((i+1):n, 1); X_io_mat_y((i+1):n, i) = X_io_mat_temp((i+1):n, 2);
    X_io_mat_x(i, (i+1):n) = X_jo_mat_temp((i+1):n, 1); X_io_mat_y(i, (i+1):n) = X_jo_mat_temp((i+1):n, 2);
    n_ij_mat_temp = n_ij{i};
    %%%%%%%%
    % n_ij_mat_temp(n - 1,:) = abs(n_ij_mat_temp(n - 1,:));
    %%%%%%%%
    n_ij_mat_x((i+1):n, i) = n_ij_mat_temp((i+1):n, 1); n_ij_mat_y((i+1):n, i) = n_ij_mat_temp((i+1):n, 2);
    n_ij_mat_x(i, (i+1):n) = -n_ij_mat_temp((i+1):n, 1); n_ij_mat_y(i, (i+1):n) = -n_ij_mat_temp((i+1):n, 2);
    [~, ind] = find(h_mat(:,i)' > 0);
    ind = ind';
    n_ij_mat = [n_ij_mat_x(:,i)'; n_ij_mat_y(:,i)'];
    F_ij = F_ij_init;%zeros(2,n);%zeros(2,length(h_ij_temp));
    T_ij = T_ij_init;%zeros(1,n);%zeros(1,length(h_ij_temp));
    %%%%%%%%
    F_grav = F_ij_init;
    F_test = F_ij_init;
    T_grav = T_ij_init;
    is_memb = ismember(ind, coll_ind_temp);
    if (sum(is_memb) > 0)
        coll_ind = unique([coll_ind, i]);
    end
    % ind = ind(ind ~= n - 1);
    is_memb = ismember(ind, coll_ind);
    n_ij_mat_per = [n_ij_mat(2,:); -n_ij_mat(1,:)];
    %%%%%%%%
    F_frot = (1-ismember(i, coll_ind))*1/2*1.2*10^(-15)*pi*(height_cell(sp_nb(i))/2).^2.*(vect_Cell_length(i) + ...
        4/6*height_cell(sp_nb(i)))*0.6*(g*t)^2.*[0;1]; %Friction force
    % F_grav(:, n - 1) = (X_io_grav(2,1) <= (height_cell(sp_nb(i))/2))*Mass_Cell(i)*g.*[0;1];
    % F_ij(:, n - 1) = (X_io_grav(2,1) <= (height_cell(sp_nb(i))/2))*(i ~= n - 1)*sqrt(0.0025/dt_forces)*4*10^4*height_cell(sp_nb(i))^1/2*(height_cell(sp_nb(i))/2 - X_io_grav(2,1))^(3/2).*[0;1];
    x = X_io_grav(:,1) - Pos_X(:,i); 
    % y = F_grav(:, n - 1);
    % T_grav(n - 1) = x(1,:).*y(2,:) - x(2,:).*y(1,:);
    T_frot = x(1,:).*F_frot(2,:) - x(2,:).*F_frot(1,:);
    % y = F_ij(:, n - 1);
    % T_ij(n - 1) = x(1,:).*y(2,:) - x(2,:).*y(1,:);
    %%%%%%%%Gravity and bottom force computation
    if ~isempty(ind)
        Bound_ind = ind > (n - nb_boundaries);
        Eff_R = height_cell(sp_nb(ind)).*height_cell(sp_nb(i))./(height_cell(sp_nb(ind)) + height_cell(sp_nb(i)));
        Eff_R(Bound_ind) = height_cell(sp_nb(i));
        %Split the force
        F_ij(:,ind) = sqrt(0.0025/dt_forces)*4*10^4*Eff_R.^(1/2)'.*h_mat(ind, i)'.^(3/2).*n_ij_mat(:,ind);     
        F_test(:,ind) = Mass_Cell(i)*is_memb'.*(n_ij_mat(2,ind) > 0).*(([0,g]*n_ij_mat(:,ind)).*n_ij_mat(:,ind)...
            + 0.8*([0,g]*n_ij_mat_per(:,ind)).*n_ij_mat_per(:,ind)...
            - 0*[0;g])/max(sum(is_memb'.*(n_ij_mat(2,ind) > 0)), 1);
        X_io_temp = [X_io_mat_x(ind, i)'; X_io_mat_y(ind, i)'];       
        x = X_io_temp - Pos_X(:,i); 
        y = F_ij(:,ind);
        y_2 = F_test(:,ind);
        T_ij(ind) = x(1,:).*y(2,:) - x(2,:).*y(1,:);
        T_grav(ind) = x(1,:).*y_2(2,:) - x(2,:).*y_2(1,:);
        ind_Pred(i) = max(mat_Pred(sp_nb(ind), sp_nb(i)))*i; %max(mat_Pred(vect_species(ind), vect_species(i)))*i; %1 if the cell should be removed due to predation, 0 otherwise
    end
    F_grav = sum(F_grav, 2) + sum(F_test, 2) - Mass_Cell(i)*g*[0;1];%(sum(is_memb'.*(n_ij_mat(2,ind) > 0)) == 0)*Mass_Cell(i)*g*[0;1];
    F_ij = sum(F_ij, 2); %2D component
    %F_ij_mex = sum_mex(F_ij,2); %2D component
    T_ij = sum(T_ij); %Radial value
    T_grav = sum(T_grav); %Radial value
    K = [fric_coeff_par*cos(cell_angle(i))^2+fric_coeff_per*sin(cell_angle(i))^2  (fric_coeff_par - fric_coeff_per)*cos(cell_angle(i))*sin(cell_angle(i)); ...
     (fric_coeff_par - fric_coeff_per)*cos(cell_angle(i))*sin(cell_angle(i)) fric_coeff_per*cos(cell_angle(i))^2+fric_coeff_par*sin(cell_angle(i))^2];
    Pos_X(:,i) = Pos_X(:,i) + K\(1./(vect_Cell_length(i) +  height_cell(sp_nb(i)) + Mass_Cell(i)).*(dt_forces*F_ij)) + 0*(F_grav + F_frot)/Mass_Cell(i)*t*dt_forces;
    cell_angle(i) = cell_angle(i) + 12/(fric_coeff_per*(vect_Cell_length(i) + 0.3*height_cell(sp_nb(i)) + Mass_Cell(i))^3)*dt_forces*T_ij + 0*0.12*dt_forces*(T_grav + T_frot)/Mass_Cell(i);
end
ind_Pred(ind_Pred == 0) = [];
%%%%%%%%Remove boundary condition
Pos_X = Pos_X(:,1:(n - nb_boundaries));
cell_angle = cell_angle(1:(n - nb_boundaries));