function props = fun_Dist_Res_Space(time_step, ind_res, rho_3D_tot, dim_Img, Name_file)
%Function to compute the distribution of the resource in space according to
%time
len_time = length(time_step);
Time = [1, round(len_time/4), round(len_time/2), 3*round(len_time/4), len_time];
[l_y, ~] = size(rho_3D_tot{ind_res, 1});
props = zeros(length(Time), l_y);
for k = 1:length(Time)
    rho_temp = rho_3D_tot{ind_res, Time(k)};
    temp = sum(rho_temp, 2);
    props(k,:) = temp/sum(temp);
    figure(k)
    bar(props(k,:)', 'stacked');
end

title_str = strcat("Resources variation according to space");
title(title_str)
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = strcat('/Users/iguex/Documents/SurfaceModels/Figures/Fig', title_str, num2str(iFig), num2str(ind_res), Name_file);
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
end
end