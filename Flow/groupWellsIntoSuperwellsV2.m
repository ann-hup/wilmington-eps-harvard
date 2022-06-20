function [superwellMatrix, xyLocationsOfWells] = groupWellsIntoSuperwellsV2(x_fluid, y_fluid, ...
    x_interval, y_interval, superwellMInput, fluidRatesFile, wellFlag, xyLocationsOfWells)
%GROUPWELLSINTOSUPERWELLSV2 Places each wilmington well into its rightful
%location in the superwellMatrix, thereby creating superwells at the end of
%this algorithm
    % Algorithm for grouping wells into superwells can be found here:
    % https://docs.google.com/document/d/1ZYKRiTQbYy7n91Eq8XhmaAIazKb1Nlvehz8pQzaIvvc/edit?usp=sharing
%  ARGUMENTS
    % wilmingtonWellIDs - list of ids corresponding with wilmington wells
    % x_fluid - list of x coordinates for wells to be grouped into superwells
    % y_fluid - list of y coordinates for wells to be grouped into superwells
    % x_interval - vector representing starting x to ending x, with a step
    % size of dx, to define the x axis partition of the grid
    % y_interval - vector representing the starting y to ending y, with a
    % step size of dy, to define the y axis partition of the grid
    % superwellMInput - inputted superwell matrix; pre-condition is that
    % the superwell matrix has been initialized
    % fluidRatesFile - matrix with fluid rates for the wells to be grouped
    % into superwells. For the valid setup of this matrix, there are n rows
    % corresponding with n time steps to the simulation, and m columns for
    % m wells to be grouped into superwells
    % wellFlag - 1 for injecting wells, -1 for producing wells - this value
    % must hold true for ALL wells being grouped here in one call of this
    % function (i.e. call this function once for injectors, and once
    % for producers)
    % xyLocationsOfWells - matrix of the UTM conversion of the well
    % locations for plotting purposes 
%  RETURNS
    % superwellMatrix - updated superwell matrix populated with wilmington
    % wells that have been consolidated into superwells, where each cell of
    % this superwellMatrix is denoted a superwell
    % xyLocationsOfWells - matrix of the UTM conversion of the well
    % locations for plotting purposes 

    
    % Group wells into their correct grid space to form superwell
    for wellPos = 1 : size(x_fluid, 1)
            % Convert WGS84 -> UTM ** USE ZONE: UTM NAD 27 **
            [UTM_x, UTM_y] = wgs2utm(y_fluid(wellPos), x_fluid(wellPos));
            
            % Get (row, col) grid location of current well
            well_col = find(UTM_x < x_interval, 1, 'first');
            % column location in superwellMatrix is wellPos_x_grid_loc - 1
            well_row_bool = flip(UTM_y < y_interval);
            well_row = find(well_row_bool == 0, 1, 'first');
            % row location in superwellMatrix is wellPos_y_grid_loc - 1
            
            well_row = (well_row - 1); well_col = (well_col - 1);

            if ~isempty(well_row) & ~isempty(well_col) & well_row > 0 & well_col > 0
                % Populate superwellMatrix accordingly
                superwellMInput(well_row, well_col).num = superwellMInput(well_row, well_col).num + 1;
                superwellMInput(well_row, well_col).wellIDs = ...
                    [superwellMInput(well_row, well_col).wellIDs; wellPos];
                xyLocationsOfWells = [xyLocationsOfWells; UTM_x, UTM_y];
                if wellFlag == 1
                    if isempty(superwellMInput(well_row, well_col).sigmaIRate)
                        superwellMInput(well_row, well_col).sigmaIRate = fluidRatesFile(:, wellPos);
                    else
                        superwellMInput(well_row, well_col).sigmaIRate = superwellMInput(well_row, well_col).sigmaIRate + ...
                            fluidRatesFile(:, wellPos);
                    end
                else
                    if isempty(superwellMInput(well_row, well_col).sigmaPRate)
                        superwellMInput(well_row, well_col).sigmaPRate = fluidRatesFile(:, wellPos);
                    else
                        superwellMInput(well_row, well_col).sigmaPRate = superwellMInput(well_row, well_col).sigmaPRate + ...
                            fluidRatesFile(:, wellPos);
                    end 
                end 
            end 
    end
    superwellMatrix = superwellMInput;
end