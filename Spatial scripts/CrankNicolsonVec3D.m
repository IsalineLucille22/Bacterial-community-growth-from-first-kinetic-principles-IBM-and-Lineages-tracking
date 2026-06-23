function [rho_fin, Mass_Res, rho_3D] = CrankNicolsonVec3D(Mass_Res, rho_3D, l_z, dy, dx, t_D, Diff_Coeff, ~, dim_Img)
dz = 1000/l_z;
l_z = l_z + 2;
[n_res, ~] = size(Mass_Res);%length(Mass_Res(:,1));
x_grid = (dim_Img(1,1) - dx/2):dx:(dim_Img(1,2) + dx/2);
y_grid = (dim_Img(2,1) - dy/2):dy:(dim_Img(2,2) + dy/2);
l_y = length(y_grid);
l_x = length(x_grid);
rho_fin = cell(1,n_res);%Put as agrugments?
rho_3D_temp = zeros(l_y, l_x, l_z);


%Concentration to change according to conditions tested and dimension of
%the simulated frame
Oxygen_Edge = 7.5342e-13;%6.5309e-14;%1.9695e-13;%9.8091e-13;%1.5695e-14;%1.1035e-11;%9.6489e-12;%9.4732e-14; %1.3491e-13; %Concentartion oxygen at the edge per resource boxe

for k = 1:n_res
    alpha = Diff_Coeff(k); %dx^2/(4*t_Diff(k)); %Thermal diffusivity
    k1 = alpha*(t_D/(dx^2));
    k2 = alpha*(t_D/(dy^2));
    k3 = alpha*(t_D/(dz^2));
    rho_3D_temp(2:(l_y-1), 2:(l_x-1), 2:(l_z-1)) = rho_3D{k};
    rho_3D_temp(2:(l_y-1), 2:(l_x-1), 2) = reshape(Mass_Res(k,:), l_y-2, l_x-2);%reshape(Mass_Res(k,:), lx-2, ly-2);%
    rho_3D_temp(rho_3D_temp<0) = 0;

%Other boundary conditions for instance for oxygen model
%     rho_mean = coeff_diff_out(k)*sum(sum(sum(rho_3D_temp(2:(lx-1), 2:(ly-1), 2:(lz-1)))))/((lx-2)*(ly-2)*(lz-2)); %Make different coefficient according to the resource %Sum the other resources and remove it in the end
%   Initial condition according to boundary conditions
    if k == 7 %&& escap_fact == 1 %If oxygen boundary condition is a reflexion, put k == 7 otherwise k == 6. 
        %Boundary conditions for oxygen model if escapment
        % rho_3D_temp(1, :, :) = Oxygen_Edge*ones(1, l_x, l_z);
        rho_3D_temp(l_y, :, :) = zeros(1,l_x,l_z);
        rho_3D_temp(:, 1, :) = rho_3D_temp(:, l_x-1, :);
        rho_3D_temp(:, l_x, :) = rho_3D_temp(:, 2, :);
        rho_3D_temp(:, :, 1) = rho_3D_temp(:, :, 2);
        rho_3D_temp(:, :, l_z) = rho_3D_temp(:, :, l_z-1);
        rho_3D_temp(1, 1, 1) = Oxygen_Edge; %Constant oxygen concentration at the l_x (lower) boundary
        rho_3D_temp(1, l_x, 1) = Oxygen_Edge;
        rho_3D_temp(l_y, l_x, 1) = 0; %No oxygen at the l_y (upper) boundary
        rho_3D_temp(l_y, 1, 1) = 0;
        rho_3D_temp(l_y, l_x, l_z) = 0;
        rho_3D_temp(1, l_x, l_z) = rho_3D_temp(2, 2, l_z-1);
        rho_3D_temp(l_y, 1, l_z) = 0;
        rho_3D_temp(1, 1, l_z) = rho_3D_temp(2, l_x-1, l_z-1);
        % Test Insulated boundary conditions
        %rho_3D_temp(l_y, :, :) = 4/3*rho_3D_temp(l_y - 1, :, :) - 1/3*rho_3D_temp(l_y - 2, :, :);
        rho_3D_temp(1, :, :) = 4/3*rho_3D_temp(2, :, :) - 1/3*rho_3D_temp(3, :, :);% Test Insulated boundary conditions 
        rho_3D_temp(1,(round(l_x/2) - round(l_x/4)):(round(l_x/2) + round(l_x/4)), 1) = Oxygen_Edge;
    else
        rho_3D_temp(1, :, :) = rho_3D_temp(l_y-1,:, :);%rho_3D_temp(2,:, :);%
        rho_3D_temp(l_y, :, :) = rho_3D_temp(2,:, :);%rho_3D_temp(l_y-1,:, :);%
        rho_3D_temp(:, 1, :) = rho_3D_temp(:, l_x-1, :);
        rho_3D_temp(:, l_x, :) = rho_3D_temp(:, 2, :);
        % Reflexive system for nutrient on x, y axis
        rho_3D_temp(:, :, 1) = rho_3D_temp(:, :, l_z-1);%rho_3D_temp(:, :, 2);%
        rho_3D_temp(:, :, l_z) = rho_3D_temp(:, :, 2);%rho_3D_temp(:, :, l_z-1);%
        rho_3D_temp(1, 1, 1) = rho_3D_temp(l_y-1, l_x-1, l_z-1);%rho_3D_temp(l_y-1, l_x-1, 2);%rho_3D_temp(2, l_x-1, 2);%
        rho_3D_temp(1, l_x, 1) = rho_3D_temp(l_y-1, 2, l_z-1);%rho_3D_temp(l_y-1, 2, 2);%rho_3D_temp(2, 2, 2);%
        rho_3D_temp(l_y, l_x, 1) = rho_3D_temp(2, 2, l_z-1);%rho_3D_temp(2, 2, 2);%rho_3D_temp(l_y-1, 2, 2);%
        rho_3D_temp(l_y, 1, 1) = rho_3D_temp(2, l_x-1, l_z-1);%rho_3D_temp(2, l_x-1, 2);%rho_3D_temp(l_y-1, l_x-1, 2);%
        rho_3D_temp(l_y, l_x, l_z) = rho_3D_temp(2, 2, 2);%rho_3D_temp(2, 2, l_z-1);%rho_3D_temp(l_y-1, 2, l_z-1);%
        rho_3D_temp(1, l_x, l_z) = rho_3D_temp(l_y-1, 2, 2);%rho_3D_temp(l_y-1, 2, l_z-1);%rho_3D_temp(2, 2, l_z-1);%
        rho_3D_temp(l_y, 1, l_z) = rho_3D_temp(2, l_x-1, 2);%rho_3D_temp(2, l_x-1, l_z-1);%rho_3D_temp(l_y-1, l_x-1, l_z-1);%
        rho_3D_temp(1, 1, l_z) = rho_3D_temp(l_y-1, l_x-1, 2);%rho_3D_temp(l_y-1, l_x-1, l_z-1);%rho_3D_temp(2, l_x-1, l_z-1);%
        % If no insulated condition for oxygen
        % rho_3D_temp(l_y, :, :) = 4/3*rho_3D_temp(l_y - 1, :, :) - 1/3*rho_3D_temp(l_y - 2, :, :);% Test Insulated boundary conditions
        % rho_3D_temp(1, :, :) = 4/3*rho_3D_temp(2, :, :) - 1/3*rho_3D_temp(3, :, :);% Test Insulated boundary conditions 
    end
    % rho_3D_temp(l_y, :, :) = 4/3*rho_3D_temp(l_y - 1, :, :) - 1/3*rho_3D_temp(l_y - 2, :, :);% Test Insulated boundary conditions
    % rho_3D_temp(1, :, :) = 4/3*rho_3D_temp(2, :, :) - 1/3*rho_3D_temp(3, :, :);% Test Insulated boundary conditions 

    rho_old = rho_3D_temp;
    rho_initial = rho_3D_temp;
    error_tol = 1;

    term1 = 1/(1+(2*k1)+(2*k2)+(2*k3));
    term2 = k1*term1;
    term3 = k2*term1;
    term4 = k3*term1;

    %starting the spatial loops.
%     rho_3D_temp(:, l_x, :) = 4/3*rho_old(:, l_x - 1, :) - 1/3*rho_old(:, l_x - 2, :);% Test Insulated boundary conditions
%     rho_3D_temp(:, 1, :) = 4/3*rho_old(:, 2, :) - 1/3*rho_old(:, 3, :);% Test Insulated boundary conditions 
    while error_tol > 1e-15
%         rho_3D_temp(2:(lx-1), 2:(ly-1), 2:(lz-1)) = rho_initial(2:(lx-1),2:(ly-1),2:(lz-1))*term1 + ...
%             (rho_old(1:(lx-2),2:(ly-1),2:(lz-1)) + rho_old(3:lx,2:(ly-1), 2:(lz-1)))*term2 + ...
%             (rho_old(2:(lx-1),1:(ly-2),2:(lz-1)) + rho_old(2:(lx-1),3:ly,2:(lz-1)))*term3 + ...
%             (rho_old(2:(lx-1),2:(ly-1),1:(lz-2)) + rho_old(2:(lx-1),2:(ly-1),3:lz))*term4; %Jacobi implicit method
        h = (rho_3D_temp(1:(l_y-2),2:(l_x-1),2:(l_z-1)) + rho_old(3:l_y, 2:(l_x-1), 2:(l_z-1)));
        v = (rho_3D_temp(2:(l_y-1),1:(l_x-2),2:(l_z-1)) + rho_old(2:(l_y-1),3:l_x,2:(l_z-1)));
        w = (rho_3D_temp(2:(l_y-1),2:(l_x-1),1:(l_z-2)) + rho_old(2:(l_y-1),2:(l_x-1),3:l_z));
        rho_3D_temp(2:(l_y-1), 2:(l_x-1), 2:(l_z-1)) = rho_initial(2:(l_y-1), 2:(l_x-1), 2:(l_z-1))*term1 + (h*term2) + (v*term3) + (w*term4); %Gauss-Seidel implicit method
        error_tol = max(max(max(abs(rho_old - rho_3D_temp))));
        rho_old = rho_3D_temp;
    end

    rho_fin{k} = rho_3D_temp(2:(l_y-1), 2:(l_x-1), 2);
    Mass_Res(k,:) = reshape(rho_fin{k}, 1, []);%reshape(rho_fin{k}', 1, []);%
    rho_3D{k} = rho_3D_temp(2:(l_y-1), 2:(l_x-1), 2:(l_z-1));
end