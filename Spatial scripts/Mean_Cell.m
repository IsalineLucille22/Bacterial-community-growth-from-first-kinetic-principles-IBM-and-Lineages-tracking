function mean_cells = Mean_Cell(cells)
[~, nb_cells] = size(cells);
mean_cells = zeros(nb_cells, 2);
for i = 1:nb_cells
    mean_cells(i, 1) = mean(cells{i});
    mean_cells(i, 2) = std(cells{i});
end
end