function vect_rate = Vect_rate(Uptake_Mat, Mu_max_Mat, Prod_Mat)
[nb_species, ~] = size(Uptake_Mat);
vect_rate = cell(nb_species, 1);
for i = 1:nb_species
    vect_rate{i} = [Uptake_Mat(i, :)' Mu_max_Mat(i, :)' Prod_Mat{i}];%Row = nb of resources, col1: kappa_1, col2: kappa_2, col3:colend: kappa_3
end
end