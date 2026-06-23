
%this script takes the final biomass or area data from simulations with the proportionality
%corrected mu to compare
% with the actual observed species biomass proportions on SE in the full minisyncom

clear
close all

species = {'Bur','Cur','Mic','Muc','Ppu','Pve','Rah'};

cd('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Tania/manuscript3/Dimalis_data/total_area_plots_trios_multi/data')

%read the file with the summed areas per species and positions

total_areas_SE = readtable('Bur-Cur-Mic-Muc-Ppu-Pve-Rah-SE-total-area.csv');

cd('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/Spatial scripts/Data');

simulations = dir('*TEST_Interappend*');%dir('*SE_mix_mu-Pos*');
% simulations = dir('*SE_mix_mu-Pos*');

close all
sz = 20;
c = orderedcolors('gem12');

%colors for the individual species here

c_species = c(find(matches(species, species)),:);

species_areas=[]; %to collect the individual summed area values per species

figH=figure;

subplot('position',[0.15,0.3,0.26,0.3]);

%add the data points from the different simulations on the same positions

for z=1:length(simulations)

    Data=load(simulations(z).name);
    
    for i=1:length(species)
    
        tmp(i)=sum(Data.Stat_Area(Data.vect_species==i));
        
        %these are colony masses
    
    end
    
    %volume of a cylinder + 2*volume of a cap = volume of a ball
    
    height_cell=Data.height_cell;
    
    length_cell=Data.length_cell;
    
    V_cell=pi*(height_cell/2).^2.*length_cell+4/3*pi*(height_cell/2).^3;
    
    % dry weights in fg of the cells
    
    m_b= 435*power(V_cell(:,1),0.86); 
    
    %dry weights in fg C of the cells
    
    m_b=m_b./2;
    
    % total masses in number of cells
    
    cells_tmp=tmp'./m_b;
    
    %cells to area
    
    cells_area=cells_tmp.*(pi*(height_cell/2).^2 + height_cell.*length_cell);
    
    for i=1:length(species)
    
        scatter(cells_area(i,1), table2array(total_areas_SE(z,i+1)),sz,c_species(i,:),'o','filled')
        hold on
    
    end
    
    
    xlim([0 1200]);
    ylim([0 1200]);
    xticks([0 200 400 600 800 1000 1200]);
    yticks([0 200 400 600 800 1000 1200]);
    
    sum_cell_area(z)=sum(cells_area(:,1));
    species_areas(z,:)=cells_area(:,1);

end
grid on

xlabel('sim sum area µm2');
ylabel('observed sum area µm2');
titlename='mu from mini syncom';
title(titlename)

% means from simulations and means from observations per species

for i=1:length(species)    
    obs(i)=mean(table2array(total_areas_SE(:,contains(total_areas_SE.Properties.VariableNames,species{i}))));
    sim(i)=mean(species_areas(:,i));
end

subplot('Position',[0.5 0.3 0.1 0.3])

stem(log2(obs./sim))
grid on
ylabel('log2 ratio obs/sim')
ylim([-2 3])
xlim([0 8])

subplot('Position',[0.70 0.3 0.1 0.3]);
bar(1:2,[mean(sum_cell_area) mean(total_areas_SE.sum)]);
hold on
scatter(1,sum_cell_area,'mo','filled')
scatter(2,total_areas_SE.sum,'ko','filled')
ylabel('sum total area µm2');
grid on


subplot('Position',[0.85 0.3 0.05 0.3]);
scatter(1,(1:7).*0.05,sz,c_species,'filled')

legendnames = species;
legend(legendnames,'Location','eastoutside');

cd('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Tania/manuscript3/Dimalis_data/Simulations')

saveas(figH,'TEST_INTER_SE-sim-obs-species-areas-mix-mu.pdf','pdf');%saveas(figH,'SE-sim-obs-species-areas-mix-mu.pdf','pdf');
