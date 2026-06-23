function [Biomass_fin, ratio, ratio_Waste, ratio_tot] = DiffBiomass_3D(Data, Res_Cons_index, Res_Prod_index)
load(Data)
Mass_Res_tot = sum(sum(sum(rho_3D_tot{Res_Cons_index, 1})));
nb_species = length(unique(vect_species));
[ratio, ratio_Waste] = deal(zeros(1, nb_species));
Biomass_fin = 0;
Biomass_init = 0;
for i = 1:nb_species
    Mass_Cell_2 = sum(Mass_Cell_Evol{i,end}) - sum(Mass_Cell_Evol{i,1});
    Biomass_fin = Biomass_fin + sum(Mass_Cell_2);
    ratio(i) = Biomass_fin/Mass_Res_tot;
    ratio_Waste(i) = sum(sum(sum(rho_3D_tot{Res_Prod_index, end})))/Mass_Res_tot;
    Biomass_init = Biomass_init + sum(Mass_Cell_Evol{i,1}); %Sum of the bacterial biomass at time 0
end
ratio_tot = (sum(sum(sum(rho_3D_tot{Res_Prod_index, end}))) + Biomass_fin + Biomass_init + sum(sum(sum(rho_3D_tot{Res_Cons_index, end}))))/(Mass_Res_tot + Biomass_init);
sum(sum(rho_3D_tot{Res_Cons_index, end}))
sum(sum(rho_3D_tot{Res_Cons_index, 1}))
