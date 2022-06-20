function nearestCellMatrix = findClosestReservoirCell(mainDir, G, superwellMatrix)
%FINDCLOSESTRESERVOIRCELL For each superwell, find the nearest cell to its
%(x, y, z) location and make sure each superwell has a unique cell location
%to penetrate when it comes to adding the wells to the W structure for the
%solver. Writes relevant data to .dat files that are to be imported for
%the W struct and schedule in addWellsAndSimulationSchedule.m
%  ARGUMENTS
    % mainDir - general directory 
    % G - meshed reservoir surface via finite element method
    % superwellMatrix - matrix populated with wilmington wells that have
    % been consolidated into superwells, where each cell of this
    % superwellMatrix is denoted a superwell
%  RETURNS
    % nearestCellMatrix - matrix containing the superwell id, its 
    % (row, col) location, chosen unique optimal cellID, and whether the
    % superwell is an injector (1) or producer (-1)
    nearestCellMatrix = [];
    takenCellIDs = [];
    superwellID = 1;
    fluid_rates_info = [];
    
    for super_row = 1 : size(superwellMatrix, 1)
        for super_col = 1 : size(superwellMatrix, 2)      
            if superwellMatrix(super_row, super_col).num > 0
                if ~isempty(superwellMatrix(super_row, super_col).sigmaIRate)
                    [nearestCellMatrix, takenCellIDs] = loopThroughCells(G, superwellMatrix, super_row, super_col, ...
                        superwellMatrix(super_row, super_col).superZ_I, 1, nearestCellMatrix, takenCellIDs, superwellID);
                    fluid_rates_info = [fluid_rates_info, superwellMatrix(super_row, super_col).sigmaIRate];
                    superwellID = superwellID + 1;
                end
                if ~isempty(superwellMatrix(super_row, super_col).sigmaPRate)
                    [nearestCellMatrix, takenCellIDs] = loopThroughCells(G, superwellMatrix, super_row, super_col, ...
                        superwellMatrix(super_row, super_col).superZ_P, -1, nearestCellMatrix, takenCellIDs, superwellID);
                    fluid_rates_info = [fluid_rates_info, superwellMatrix(super_row, super_col).sigmaPRate];
                    superwellID = superwellID + 1;
                end
            end
        end
    end
    
    % Write superwell id, chosen unique optimal cellID, and well type
    superwell_cell_info = nearestCellMatrix(:, [1, 4, 5]);
    well_cell_info_file = fopen(strcat(mainDir,'superwellData/well_cell_mapping.dat'),'w');
    fprintf(well_cell_info_file, '%d %d %d\n', superwell_cell_info.');
    
    % Convert fluid rates from barrels/month to m^3/s
    % 1 bbl = 0.158987 m^3
    % 1 month = 2.628e+6 s
    % OR 1 fluid barrel/month = 4.53429987e-8 m^3/s
    fluid_rates_info = fluid_rates_info;
    fluid_rates_info = fluid_rates_info * 4.53429987e-8;
    
    % Write fluid rates, where each column corresponds with a superwell and
    % its column vector is the fluid rates, with one rate for each timestep
    writematrix(fluid_rates_info, strcat(mainDir,'superwellData/superwell_fluid_rates.dat'));
end


function [nearestCellMatrix, takenCellIDs] = loopThroughCells(G, superwellMatrix, ...
    super_row, super_col, z_depth, wellSign, nearestCellMatrix, takenCellIDs, superwellID)
    bestDistCell = -1; bestCellID = -1;
    for cellID = 1 : length(G.cells.centroids)
        if ~ismember(G.cells.indexMap(cellID), takenCellIDs)
            cur_cell_centr = G.cells.centroids(cellID, :);
            dist_to_cell = sqrt((superwellMatrix(super_row, super_col).superX - cur_cell_centr(1))^2 + ...
                            (superwellMatrix(super_row, super_col).superY - cur_cell_centr(2))^2 + ...
                            (z_depth - cur_cell_centr(3))^2);

            if dist_to_cell < bestDistCell | bestDistCell == -1
                bestDistCell = dist_to_cell; bestCellID = G.cells.indexMap(cellID);
            end
        end
    end
    
    nearestCellEntry = [superwellID, super_row, super_col, bestCellID, wellSign];
    nearestCellMatrix = [nearestCellMatrix; nearestCellEntry];
    takenCellIDs = [takenCellIDs; bestCellID];   
end