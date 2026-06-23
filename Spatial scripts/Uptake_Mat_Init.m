function Uptake_Mat = Uptake_Mat_Init(mu_max, kappa1_mat, kappa3_mat, nb_species, nb_resources)
% Default initialization of the uptake matrix assumes that each species
% consumes one resource (the first one).
% Otherwise, manually initialize the production matrix.
Uptake_Mat = zeros(nb_species, nb_resources);
for i = 1:nb_species
    Uptake_Mat(i, 1) = (mu_max(i) + kappa3_mat(i))/kappa1_mat(i);
end
end