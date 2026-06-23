function Prod_Mat = Init_Prod_Mat(kappa_mat, nb_species, nb_resources, pos_main_waste)
% Default initialization of the production matrix assumes that each species
% consumes one resource and produces the same unique waste product.
% Otherwise, manually initialize the production matrix.
Prod_Mat = cell(nb_species, 1);
for i = 1:nb_species
    Prod_Mat{i} = zeros(nb_resources);
    Prod_Mat{i}(1, pos_main_waste) = kappa_mat(i);
end
end