function cellNeighbors = getCellNeighbors(G, cellID) 
    cell_n_locs = any(G.faces.neighbors == cellID, 2);
    cell_n_matrix = G.faces.neighbors(cell_n_locs, :);
    cellNeighbors = [];
    for entry = 1 : numel(cell_n_matrix) 
        if cell_n_matrix(entry) ~= 0 && cell_n_matrix(entry) ~= cellID
            cellNeighbors = [cellNeighbors; cell_n_matrix(entry)];
        end
    end
end
