function Mat_rate = Adjust_Mat_rate(kappa_3, Mat_rate)
%This function allows to adjust the Mat_rate matrices per cells according
%to the factor multiplication applied to kappa_3. This will proceed per
%cell and per resource. This function is useful only if we want to modify
%individual cell yields, for instance due to spatial interactions.
[nb_res, ~] = size(kappa_3);
for i = 1:nb_res
    temp_kappa_3 = kappa_3(i,:);
    temp_Mat_rate = squeeze(Mat_rate(i, 3:end, :));
    fact_mod = temp_kappa_3./sum(temp_Mat_rate, 1);
    Mat_rate(i, 3:end, :) = Mat_rate(i, 3:end, :).*reshape(fact_mod,1, 1, []);
    Mat_rate(isinf(Mat_rate) | isnan(Mat_rate)) = 0;
end
end