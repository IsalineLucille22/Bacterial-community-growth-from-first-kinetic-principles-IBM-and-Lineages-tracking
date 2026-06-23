function [rho_fin, Mass_Res, rho_3D] = BoltzmannLattice(Mass_Res, dx, t_D, t_Diff, t_LB, id_tot, dim_Img)
[n_res, ~] = size(Mass_Res);% length(Mass_Res(:,1));
x_grid = (dim_Img(1,1) + dx/2):dx:(dim_Img(1,2));
y_grid = (dim_Img(2,1) + dx/2):dx:(dim_Img(2,2)); %Change dy
lx = length(x_grid);
ly = length(y_grid);
n_step = round(t_D/t_LB);
% We put 0 for the weight corresponding to d = 0, but only to have zeros in
% the first column of fOut, for computation. We should have the sum of
% weight equals to 1, so we have w_0 = 4/9 according to this.
w  = [0; 1/4; 1/4; 1/4; 1/4]; %[0; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36]; %[4/9; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36]; %Weight vector
scaling_time = t_LB./t_Diff;

non_zero_res = find(sum(Mass_Res,2) > 0);
rho_fin = cell(1,n_res);

for i = 1:n_step
    rho = reshape(Mass_Res, n_res, lx, ly);
    rho_fin = cell(1,n_res);
%     non_zero_res = find(sum(Mass_Res,2) > 0);
    
    
    
    % COLLISION STEP %Particles in the same site interact  
    for k = 1:length(non_zero_res)
        ind_temp = non_zero_res(k);
        scaling_time_temp = scaling_time(ind_temp);
        fOut = scaling_time_temp*w.*rho(ind_temp,:,:);%.*Proba_free;%
        fOut(1,:,:) = rho(ind_temp,:,:) - sum(fOut);
        
        fOut = fOut(id_tot);
        
        % MACROSCOPIC VARIABLES
        rho_temp = sum(fOut); %sum(fIn); %Sum on the 9 directions
%         rho_temp(rho_temp<0) = 0;
        temp = reshape(rho_temp, lx, ly);%squeeze(rho_temp);
        rho_fin{ind_temp} = temp;%reshape(rho_temp, lx, ly);%
%         rho_3D{ind_temp} = rho_fin{ind_temp};
        Mass_Res(ind_temp,:) = temp(:)';%reshape(rho_fin{ind_temp}, 1, []);
    end
end
rho_3D = rho_fin;