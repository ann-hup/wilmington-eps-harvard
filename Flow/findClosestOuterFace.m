function [superwellMatrix, bestFaceIDs] = findClosestOuterFace(superwellMInput, upperExteriorFaceCentroids, ...
    upperExteriorFaceIDs)
%FINDCLOSESTOUTERFACE For each superwell, find the closest upper exterior
%boundary face, and then use the ultimately chosen optimal face id to set
%the preliminary z depth for that superwell. The final z depth at this
%stage in the algorithm that is set for the superwell is the subtraction of
%a fixed depth depending on whether the well is an injector or producer. 
    % ARGUMENTS
    % superwellMInput - inputted superwell matrix; pre-condition is that the
    % wilmington wells have already been grouped into their corresponding
    % superwell
    % upperExteriorFaceCentroids - c x 3 matrix, where c is the number of
    % upper exterior face ids found (length(upperExteriorFaceIDs) = c),
    % with the centroid (x, y, z) of each face
    % upperExteriorFaceIDs - list of ids corresponding with upper exterior
    % boundary faces on the top of the rock reservoir
    % RETURNS
    % superwellMatrix - updated superwell matrix; post-condition is that
    % the initial z depth of each superwell has been set
    % bestFaceIDs - list of chosen face ids in this algorithm that were
    % picked by each superwell as its closest 2D projection on top of the
    % rock
    bestFaceIDs = [];
    % For each superwell, find it's closest 2D projection by finding its
    % closest face on the upper exterior boundary of the rock reservoir
    for super_row = 1 : size(superwellMInput, 1)
        for super_col = 1 : size(superwellMInput, 2)          
            if superwellMInput(super_row, super_col).num > 0
                bestDist = -1; bestFace = -1; bestFaceID = -1;
                for upperFaceCentr = 1 : length(upperExteriorFaceCentroids)
                    cur_centr = upperExteriorFaceCentroids(upperFaceCentr, :);
                    dist_to_face = sqrt((superwellMInput(super_row, super_col).superX - cur_centr(1))^2 + ...
                        (superwellMInput(super_row, super_col).superY - cur_centr(2))^2);
                    if dist_to_face < bestDist | bestDist == -1
                        bestDist = dist_to_face; bestFace = cur_centr; 
                        bestFaceID = upperExteriorFaceIDs(upperFaceCentr);
                    end
                end
                
                bestFaceIDs = [bestFaceIDs; bestFaceID];
                
                % Set the depth of the superwell to be a fixed distance
                % from closest face centroid on upper exterior boundary.
                % Injectors are believed to be displaced at a greater depth
                % than producers
                if ~isempty(superwellMInput(super_row, super_col).sigmaIRate)
                    superwellMInput(super_row, super_col).superZ_I = bestFace(3) - 800;
                end
                if ~isempty(superwellMInput(super_row, super_col).sigmaPRate)
                    superwellMInput(super_row, super_col).superZ_P = bestFace(3) - 200;
                end                
            end
        end
    end
    superwellMatrix = superwellMInput;
end