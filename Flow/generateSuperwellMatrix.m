function [wilmingtonWellIDs, upperExteriorFaceIDs, upperExteriorFaceCentroids, superwellMatrix] = ...
    generateSuperwellMatrix(mainDir, G, x_I, y_I, x_P, y_P, injectionRates, productionRates)
%GENERATESUPERWELLMATRIX Generates the superwell matrix by grouping every possible well into
%the cell it belongs to, based on the stratification grid that is
%devised by bounds for x and y and a chosen dx and dy. Each cell is then
%treated as a superwell (note that each cell may have both an injector
%superwell and a producer superwell). Each superwell has a (x, y, z)
%location and a vector representing the superwell's injection or 
%production fluid rates at each time step, to capture the total behavior of
%all the wells that were assigned its particular cell.
    % - Get ids of wells that are wilmington wells
    % - Specify the partition grid properties 
    % - Build the empty superwell matrix
    % - Group wells into superwells according to an algorithm devised in a
    % separate function groupWellsIntoSuperwells. Algorithm theory
    % can be found here: 
    % https://docs.google.com/document/d/1ZYKRiTQbYy7n91Eq8XhmaAIazKb1Nlvehz8pQzaIvvc/edit?usp=sharing
%  ARGUMENTS
    % mainDir - general directory 
    % G - meshed reservoir surface via finite element method
    % x_I - list of x coordinates for all injecting wells in original .dat file
    % y_I - list of y coordinates for all injecting wells in original .dat file
    % x_P - list of x coordinates for all producing wells in original .dat file
    % y_P - list of y coordinates for all producing wells in original .dat file
    % injectionRates - matrix with injection rates, where there are n rows
    % corresponding with n time steps to the simulation, and m columns for
    % m wells in the original .dat file
    % productionRates - matrix with production rates, where there are n rows
    % corresponding with n time steps to the simulation, and m columns for
    % m wells in the original .dat file
%  RETURNS
    % wilmingtonWellIDs - list of ids corresponding with wilmington wells
    % upperExteriorFaceIDs - list of ids corresponding with upper exterior
    % boundary faces on the top of the rock reservoir
    % upperExteriorFaceCentroids - c x 3 matrix, where c is the number of
    % upper exterior face ids found (length(upperExteriorFaceIDs) = c),
    % with the centroid (x, y, z) of each face
    % superwellMatrix - matrix populated with wilmington wells that have
    % been consolidated into superwells, where each cell of this
    % superwellMatrix is denoted a superwell

    % Get wilmington well IDs
    wilmingtonWellIDs = getWilmingtonWells(mainDir);
    
    % Get the ids and centroids of the upper exterior faces on the
    % reservoir, as these will come in handy when finding the closest face
    % for each superwell on the top of the rock
    [upperExteriorFaceIDs, upperExteriorFaceCentroids] = getExtremaCellsAndFaces(G);
    
    % Define partition grid specifications to group wells into superwells
    
    % Partition properties when z axis is reversed (0 - top, < 0 - bottom)
    % dx = 0.01e05; dy = 0.001e06;
    % start_x = 3.8e05; end_x = 3.95e05; % 15
    % start_y = 3.732e06; end_y = 3.742e06; % 10 --> 150 possible superwells
    
    % Partition properties when z axis is normal (0 - bottom, > 0 - top)
     dx = 0.01e05 / 5; dy = 0.001e06 / 5;
     start_x = 3.76e05; end_x = 4.01e05;
     start_y = 3.727e06; end_y = 3.752e06; 
    
    x_interval = start_x : dx: end_x; y_interval = start_y : dy : end_y;
    reverse_y_interval = flip(y_interval);
   
    % Initalize the empty superwell matrix, where each cell is a superwell
    superwellMatrix = initializeSuperwellMatrix(x_interval, y_interval, reverse_y_interval);
    
    % Place each wilmington well into its rightful cell, thereby
    % consolidating that well into what will eventually be its
    % corresponding superwell
    %[superwellMatrix, xyLocationsOfWells] = groupWellsIntoSuperwells(wilmingtonWellIDs, x_I, y_I, x_interval, y_interval, ...
                                              %  superwellMatrix, injectionRates, productionRates);
    [superwellMatrix, xyLocationsOfWells] = groupWellsIntoSuperwellsV2(x_I, y_I, x_interval, y_interval, ...
                                                superwellMatrix, injectionRates, 1, []);
    [superwellMatrix, xyLocationsOfWells] = groupWellsIntoSuperwellsV2(x_P, y_P, x_interval, y_interval, ...
                                                superwellMatrix, productionRates, -1, xyLocationsOfWells);                                        
    makeSuperwellGrid(G, xyLocationsOfWells);
end


function wilmingtonWells = getWilmingtonWells(mainDir)
%GETWILMINGTONWELLS Get the corresponding well ids of wilmington wells
    fieldNameFile = fopen(strcat(mainDir, 'dataExported/Injection/fieldName.dat'));
    fieldNames = textscan(fieldNameFile, '%s', 'delimiter', '\n');
    fclose(fieldNameFile);
    fieldNames = string(fieldNames{:});
    wilmingtonWells = find(fieldNames == "Old Wilmington (ABD)" | fieldNames == "Wilmington");
end


function [upperExteriorFaceIDs, upperExteriorFaces] = getExtremaCellsAndFaces(G)
%GETEXTREMACELLSANDFACES Get the face ids and corresponding centroids of
%the upper exterior boundary faces, which will eventually become
%candidates for each superwell's closest face id on the top of the rock
    exteriorFaceLocs = any(G.faces.neighbors == 0, 2); % bool vector
    exteriorFaceCentroids = G.faces.centroids(exteriorFaceLocs, :);
    exteriorFaceIDs = find(exteriorFaceLocs == 1);
    
    % Get upper exterior boundary cells and faces only with this cutoff
    % depth to filter properly
    cutoffDepthFace = abs((min(G.faces.centroids(:, 3)) + max(G.faces.centroids(:, 3))) / 2);
    % Note: 3799 cells in this subset
    
    valid_positions = any(abs(exteriorFaceCentroids(:, 3)) < 1200, 2);
    upperExteriorFaces = exteriorFaceCentroids(valid_positions, :);
    upperExteriorFaceIDs = exteriorFaceIDs(valid_positions);
    % Note: 3921 faces in this subset
    % Note: 2161 faces in the subset with cutoff depth fixed at abs(-1200)
end


function superwellMatrix = initializeSuperwellMatrix(x_interval, y_interval, reverse_y_interval)
%INITIALIZESUPERWELLMATRIX Initalize the superwell matrix, where each cell
%is to be treated as a superwell. Each cell has the properties noted in the
%for loop below.
    superwellMatrix(1 : (length(y_interval) - 1), 1 : (length(x_interval) - 1)) = struct;
    for r = 1 : (length(y_interval) - 1)
        for c = 1 : (length(x_interval) - 1)
            superwellMatrix(r, c).num = 0;
            superwellMatrix(r, c).wellIDs = [];
            superwellMatrix(r, c).sigmaIRate = [];
            superwellMatrix(r, c).sigmaPRate = [];
            superwellMatrix(r, c).superX = (x_interval(c) + x_interval(c + 1)) / 2;
            superwellMatrix(r, c).superY = (reverse_y_interval(r) + reverse_y_interval(r + 1)) / 2;
            superwellMatrix(r, c).superZ_I = 0;
            superwellMatrix(r, c).superZ_P = 0;
        end
    end
end