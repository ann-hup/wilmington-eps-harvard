function [schedule, W] = Create_Well_and_SimulationSchedule(mainDir, G, bc, rock)
%% Function summary

    % Find the corresponding cellID of the optimal cell for the placement of each well
    
    % Only filter for Old Wilmington ABD or Wilmington
    fieldNameFile = fopen(strcat(mainDir, 'dataExported/Injection/fieldName.dat'));
    fieldNames = textscan(fieldNameFile, '%s', 'delimiter', '\n');
    fclose(fieldNameFile);
    fieldNames = string(fieldNames{:});
    wilmingtonWellIDs = find(fieldNames == "Old Wilmington (ABD)" | fieldNames == "Wilmington");
    % Note: 7316 wells in this subset
    
    % Get cells and faces at exterior boundaries (upper and lower)
    exteriorFaceLocs = any(G.faces.neighbors == 0, 2); % bool vector
    exteriorFaces = G.faces.neighbors(exteriorFaceLocs, :);
    exteriorFacesVec = reshape(exteriorFaces, [numel(exteriorFaces), 1]);
    exteriorCellIDs = unique(exteriorFacesVec(exteriorFacesVec ~= 0));
    exteriorFaceCentroids = G.faces.centroids(exteriorFaceLocs, :);
    
    % Get upper exterior boundary cells and faces only (using cutoff depths as
    % filter for this demarcation)
    cutoffDepthCell = (min(G.cells.centroids(:, 3)) + max(G.cells.centroids(:, 3))) / 2;
    cutoffDepthFace = (min(G.faces.centroids(:, 3)) + max(G.faces.centroids(:, 3))) / 2;
    upperExteriorCellIDs =  exteriorCellIDs(arrayfun(@(cellID) G.cells.centroids(cellID, 3) > cutoffDepthCell, exteriorCellIDs));
    % Note: 3799 cells in this subset
    
    upperExteriorFaceCentroids = exteriorFaceCentroids(any(exteriorFaceCentroids(:, 3) > cutoffDepthFace, 2), :);
    % Note: 3921 faces in this subset
        
    
    % Overall objective: for each superwell, find nearest cell from upper exterior boundary
    
    % Read in coordinate data for injectors and producers
    injectionWellsFile = dlmread(strcat(mainDir, 'dataExported/Injection/coord.dat'));  
    productionWellsFile = dlmread(strcat(mainDir, 'dataExported/Production/coord.dat'));  
    x_I = injectionWellsFile(:, 2); y_I = injectionWellsFile(:, 1);
    x_P = productionWellsFile(:, 2); y_P = productionWellsFile(:, 1);
    % Note: right now, x_I = x_P and y_I = y_P
    
    % Read in injection/production rates
    injectionRatesFile = dlmread(strcat(mainDir, 'dataExported/Injection/waterInjection.dat'));  
    productionRatesFile = dlmread(strcat(mainDir, 'dataExported/Production/waterProduction.dat')); 
    
    
    % First, need to create superwells by iterating over wells according to
    % grid created with dx and dy (see makeSuperwells.m) and grouping them
    
    % Specifications of partitioned grid to group wells into superwells
    dx = 0.01e05; dy = 0.001e06;
    start_x = 3.8e05; end_x = 3.95e05;
    start_y = 3.732e06; end_y = 3.742e06;
    x_interval = start_x : dx: end_x; y_interval = start_y : dy : end_y;
    reverse_y_interval = flip(y_interval);
    
    % Initialize tracker matrices
    validWellIDs = []; bestCellIDs = []; newData = []; invalidWellIDs = [];
    nonzeroRateWellIDs = [];
    % valid wells has 4636 wells, but ends up being 4632 with 4 discounted
    
    % Initialize superwell matrix
    superwellMatrix(1 : (length(y_interval) - 1), 1 : (length(x_interval) - 1)) = struct;
    for iter1 = 1 : (length(y_interval) - 1)
        for iter2 = 1 : (length(x_interval) - 1)
            superwellMatrix(iter1, iter2).num = 0;
            superwellMatrix(iter1, iter2).wellIDs = [];
            superwellMatrix(iter1, iter2).sigmaIRate = [];
            superwellMatrix(iter1, iter2).sigmaPRate = [];
            superwellMatrix(iter1, iter2).superX = (x_interval(iter2) + x_interval(iter2 + 1)) / 2;
            superwellMatrix(iter1, iter2).superY = (reverse_y_interval(iter1) + reverse_y_interval(iter1 + 1)) / 2;
            superwellMatrix(iter1, iter2).superZ_I = 0;
            superwellMatrix(iter1, iter2).superZ_P = 0;
            % may potentially add more fields here ...
        end
    end
    
    % can get rid of this construct later 
    tempMatrix = [];
    
    % Group wells into their correct grid space to form superwell
    for wellPos = 1 : length(wilmingtonWellIDs)
        rateSum_I = sum(injectionRatesFile(:, wilmingtonWellIDs(wellPos)), 1);
        rateSum_P = sum(productionRatesFile(:, wilmingtonWellIDs(wellPos)), 1);
        if (rateSum_I == 0 & rateSum_P == 0) | (rateSum_I ~= 0 & rateSum_P ~= 0)
            % means well is invalid
            invalidWellIDs = [invalidWellIDs; wilmingtonWellIDs(wellPos)];
            
            % see how many wells have nonzero rates for both I/P
            if (rateSum_I ~=0 & rateSum_P ~= 0)
                nonzeroRateWellIDs = [nonzeroRateWellIDs; wilmingtonWellIDs(wellPos)];
            end
            % need to consider if this invalid well conditional is
            % necessary ... because adding 0 vector doesn't make sigmaIRate
            % or sigmaPRate incorrect ... but are there wells with both
            % nonzero rates?
        else
            % means wells is valid
            validWellIDs = [validWellIDs; wilmingtonWellIDs(wellPos)];
            
            % Convert WGS84 -> UTM
            [UTM_x, UTM_y] = wgs2utm(y_I(wilmingtonWellIDs(wellPos)), x_I(wilmingtonWellIDs(wellPos)));
            
            % Get (row, col) grid location of current well
            well_col = find(UTM_x < x_interval, 1, 'first');
            % column location in superwellMatrix is wellPos_x_grid_loc - 1
            well_row_bool = flip(UTM_y < y_interval);
            well_row = find(well_row_bool == 0, 1, 'first');
            % row location in superwellMatrix is wellPos_y_grid_loc - 1
            
            % Handle wells not in domain of (x_interval) x (y_interval)
            well_row = (well_row - 1); well_col = (well_col - 1);
            %{
            if isempty(well_row) | well_row <= 0
                well_row = length(y_interval) - 1;
            end
            if isempty(well_col) | well_col <= 0
                well_col = 1;
            end
            %}
            if ~isempty(well_row) & ~isempty(well_col) & well_row > 0 & well_col > 0
                tempEntry = [wilmingtonWellIDs(wellPos), well_row, well_col];
                tempMatrix = [tempMatrix; tempEntry];
                % we should eliminate the four wells not in the main grid


                % Populate superwellMatrix accordingly
                superwellMatrix(well_row, well_col).num = superwellMatrix(well_row, well_col).num + 1;
                superwellMatrix(well_row, well_col).wellIDs = ...
                    [superwellMatrix(well_row, well_col).wellIDs; wilmingtonWellIDs(wellPos)];
                if rateSum_I ~= 0
                    if isempty(superwellMatrix(well_row, well_col).sigmaIRate)
                        superwellMatrix(well_row, well_col).sigmaIRate = injectionRatesFile(:, wilmingtonWellIDs(wellPos));
                    else
                        superwellMatrix(well_row, well_col).sigmaIRate = superwellMatrix(well_row, well_col).sigmaIRate + ...
                            injectionRatesFile(:, wilmingtonWellIDs(wellPos));
                    end

                else
                    if isempty(superwellMatrix(well_row, well_col).sigmaPRate)
                        superwellMatrix(well_row, well_col).sigmaPRate = productionRatesFile(:, wilmingtonWellIDs(wellPos));
                    else
                        superwellMatrix(well_row, well_col).sigmaPRate = superwellMatrix(well_row, well_col).sigmaPRate + ...
                            productionRatesFile(:, wilmingtonWellIDs(wellPos));
                    end 
                end  
            end
        end
    end
    
    % For each superwell, find it's closest 2D projection by finding its
    % closest face on the upper exterior boundary of the rock reservoir
    for super_row = 1 : size(superwellMatrix, 1)
        for super_col = 1 : size(superwellMatrix, 2)          
            % superWell = superwellMatrix(super_row, super_col);
            if superwellMatrix(super_row, super_col).num > 0
                bestDist = -1; bestFace = -1;
                for upperFaceCentr = 1 : length(upperExteriorFaceCentroids)
                    cur_centr = upperExteriorFaceCentroids(upperFaceCentr, :);
                    dist_to_face = sqrt((superwellMatrix(super_row, super_col).superX - cur_centr(1))^2 + ...
                        (superwellMatrix(super_row, super_col).superY - cur_centr(2))^2);
                    if dist_to_face < bestDist | bestDist == -1
                        bestDist = dist_to_face; bestFace = cur_centr;
                    end
                end
                % Set the depth of the superwell to be a fixed distance
                % from closest face centroid on upper exterior boundary.
                % Injectors are displaced at a greater depth than
                % producers
                if ~isempty(superwellMatrix(super_row, super_col).sigmaIRate)
                    superwellMatrix(super_row, super_col).superZ_I = bestFace(3) - 800;
                end
                if ~isempty(superwellMatrix(super_row, super_col).sigmaPRate)
                    superwellMatrix(super_row, super_col).superZ_P = bestFace(3) - 200;
                end                
            end
        end
    end
    
    % For each superwell, find the nearest cell to its (x, y, z) location
    % and make sure each superwell has a unique cell location to penetrate
    nearestCellMatrix = [];
    
    nearestCellMatrix = findClosestReservoirCell(mainDir, G, superwellMatrix);
    
    for super_row = 1 : size(superwellMatrix, 1)
        for super_col = 1 : size(superwellMatrix, 2)      
            superWell = superwellMatrix(super_row, super_col);
            if superwellMatrix(super_row, super_col).num > 0
                if ~isempty(superwellMatrix(super_row, super_col).sigmaIRate)
                    bestDistCell_I = -1; bestCellID_I = -1;
                    for cellID_I = 1 : length(G.cells.centroids)
                        cur_cell_centr = G.cells.centroids(cellID_I, :);
                        dist_to_cell_I = sqrt((superwellMatrix(super_row, super_col).superX - cur_cell_centr(1))^2 + ...
                        (superwellMatrix(super_row, super_col).superY - cur_cell_centr(2))^2 + ...
                        (superwellMatrix(super_row, super_col).superZ_I - cur_cell_centr(3))^2);
                        if dist_to_cell_I < bestDistCell_I | bestDistCell_I == -1
                            bestDistCell_I = dist_to_cell_I; bestCellID_I = G.cells.indexMap(cellID_I);
                        end
                    end
                    nearestCellEntry = [super_row, super_col, bestCellID_I, 1];
                    nearestCellMatrix = [nearestCellMatrix; nearestCellEntry];
                end
                if ~isempty(superwellMatrix(super_row, super_col).sigmaPRate)
                    bestDistCell_P = -1; bestCellID_P = -1;
                    for cellID_P = 1 : length(G.cells.centroids)
                        cur_cell_centr = G.cells.centroids(cellID_P, :);
                        dist_to_cell_P = sqrt((superwellMatrix(super_row, super_col).superX - cur_cell_centr(1))^2 + ...
                        (superwellMatrix(super_row, super_col).superY - cur_cell_centr(2))^2 + ...
                        (superwellMatrix(super_row, super_col).superZ_P - cur_cell_centr(3))^2);
                        if dist_to_cell_P < bestDistCell_P | bestDistCell_P == -1
                            bestDistCell_P = dist_to_cell_P; bestCellID_P = G.cells.indexMap(cellID_P);
                        end
                    end 
                    nearestCellEntry = [super_row, super_col, bestCellID_P, -1];
                    nearestCellMatrix = [nearestCellMatrix; nearestCellEntry];
                end
            end
        end
    end
                                    % --TO DO -- % 
                                    
    % NEXT STEP: MAKE CELL LOCATIONS UNIQUE
    % GO THROUGH SUPERWELLS AND ADD TO W STRUCT
    % FIX CODE - MODULARIZE/ABSTRACT INTO DIFFERENT FUNCTION FILES
    %   MAKE SURE TO SAVE SUPERWELL DATA TO .DAT FILE CONTAINING (X, Y, Z)
    %   LOCATION OF SUPERWELL AND ANY OTHER RELEVANT DATA (I.E.
    %   PRODUCER/INJECTOR, ANY CELL/FACE/WELL IDS, ETC.)
    
    % nearestCellMatrix is 96x4, so there are 96 superwells
    % 47 are injectors, 49 are producers
    % 84 cell locations are unique, so there are 12 wells that are in the
    % location as another superwell --> how to fix this? 
    % perhaps a visited list 
    % if cell has been visited, find the next nearest cell and see if that
    % has been visited already
 
    % next task: add super wells to struct W
    % go through simple for loop of super wells matrix
    
    % next task: create scheduler
    % 1 month, runs from 1938-2019 - currently just 1978-2019
    % ask more about this...
    
     % FIGURE OUT EVAL CACHE ISSUE
            
            % THEN USE THIS MATRIX WITH FOR LOOP CODE BELOW TO FIND NEAREST
            % CELL --> CREATE NEEDED DEPTH --> ADD SUPERWELLS TO W STRUCT
            % --> CREATE SCHEDULE
            
                               % -- TO DO -- %
   
                             % -- OLDER CODE -- % 
    %{  
    SOME REMNANT OF ADDWELL CODE
        best_dist_LIST(w_index) = best_dist;
        best_cell_index_LIST(w_index) = G.cells.indexMap(best_cell_index, :);
        W = addWell(W,G,rock, best_cell_index, 'type','rate', 'val', ...
            injectionRate(1, w_index), 'Comp_i', [1 0], 'Sign', well_type(w_index), 'refDepth', -3000 );
    
    % Write optimal CellIDs to file
    wellCellIDFile = fopen(strcat(mainDir,'fromJulia/WellCellId.dat'),'w');
    fprintf(wellCellIDFile, '%d\n', best_cell_index_LIST.');
    %}

    %% Computing the injection rate here
    %fractionOfPoreVolumeToInject = 0.02;
    %injectionTime = 50*365*86400; %20 years;
    %PV = sum(G.cells.volumes .* G.cells.cellPoro);  %% Total pore volume space
    %injectedCO2PoreVolume = PV * fractionOfPoreVolumeToInject
    %injectionRate = injectedCO2PoreVolume  / injectionTime
    %injectionRate = 1*17.6 % 1 Mt / year assuming a surface density of 1.8
    %injectionRate  = 1e-6
    %massOfCO2Injected_MT = injectedCO2PoreVolume * 1.8 / 1e9

    % tmp = dlmread('fromJulia/InjectionRatesToMRST.dat');
    % dt = tmp(:,1);
    % injectionRate = tmp(:,2);
    %% Create wells
    % W = [];
    % %W = addWell(W,G,rock, wellCellId, 'type','rate', 'val', 0.00000385, 'Comp_i', [1 0] );
    %for k = 1:length(dt)
        %W = addWell(W,G,rock, wellCellId, 'type','rate', 'val', injectionRate(k) , 'Comp_i', [0 1] );
    %end
    
    % W = addWell(W,G,rock, wellCellId, 'type','rate', 'val', injectionRate(1) , 'Comp_i', [1 0], 'Sign', 1, 'refDepth', -3000 );

    W = [];
    %% Create schedule here
    factorWI = 0.2
    %factorWI = 1
    schedule=[];
 
    schedule = simpleSchedule(dt,'W',W(1), 'bc', bc);
    schedule.step.control=[1:length(dt)]';
    
    scheduleJS=[];
    tmp = cell(length(dt), 1);  %Create schedules
    scheduleJS = struct('step', schedule.step);
    scheduleJS.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
    for wellID = 1:length(W)
        for k=1:length(dt)
            scheduleJS.control(k).W(wellID) = W(wellID);
            scheduleJS.control(k).W(wellID).val = injectionRate(k, wellID);
            %scheduleJS.control(k).W.lims.bhp = 40e6;
            %scheduleJS.control(k).W.WI = scheduleJS.control(k).W.WI*factorWI;
        end
    end
    
    schedule=scheduleJS;
end