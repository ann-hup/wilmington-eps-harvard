function [schedule, W] = TwoPhaseADWells(mainDir, G, bc, rock)
%TWOPHASEADWELLS Summary of this function goes here
%   Detailed explanation goes here
% CODE HERE

% PLAN:
% plot wells + G struct to pick good wells to test

% plot initial saturation (follow pressure plot mechanism)

% modify in this function with 1 injector/1 producer to simulate schedule
% with --> then run simulation and see how solver reacts with this base
% case (do all timesteps run? anything run over maximum iterations?)

    well_cell_data = dlmread(strcat(mainDir, 'superwellData/well_cell_mapping.dat'));  
    well_fluid_rate_data = dlmread(strcat(mainDir, 'superwellData/superwell_fluid_rates.dat')); 
    
    % import timestep data
   % timestep_info = dlmread(strcat(mainDir, 'dataExported/Injection/date.dat'));  
  %  dt = diff(timestep_info);
  %  dt = dt * 3.154e7;
    dt = repelem(30 * day, 20);
    
     W = [];
     % maxX position from Well cells is cell 3658, located at Well 1655 -
     % injector
     
     % minX position form Well cells is cell 17168, located at Well 11 -
     % producer
     
     
    % W = addWell(W, G, rock, 3658, 'type', 'rate', 'val', ...
    %     well_fluid_rate_data(1, 1655), 'Comp_i', [1 0], 'Sign', 1, 'refDepth', []);
     
     W = addWell(W, G, rock, 17168, 'type', 'rate', 'val', ...
         -1 * well_fluid_rate_data(1, 11), 'Comp_i', [0 1], 'Sign', -1, 'refDepth', 0);
     
   %  W = addWell(W, G, rock, well_cell_data(entry, 2), 'type','rate', 'val', ...
    %          well_fluid_rate_data(1, entry), 'Comp_i', [1 0], 'Sign', 1, 'refDepth', []); % [0, 1]
    
    schedule = [];
    factorWI = 0.2;
     
    schedule = simpleSchedule(dt,'W',W, 'bc', bc);
    schedule.step.control=[1:length(dt)]';
    
    scheduleJS = [];
    
    tmp = cell(length(dt), 1); 
    scheduleJS = struct('step', schedule.step);
    scheduleJS.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
    
   % W_WellIDs = [1655, 11];
   W_WellIDs = [11];
   
   n_controls = 20;
    
    for wellID = 1:length(W)
        for k = 1:n_controls %length(dt)
            scheduleJS.control(k).W(wellID) = W(wellID);        
            if W(wellID).sign == 1
                % scheduleJS.control(k).W(wellID).val = well_fluid_rate_data(k, W_WellIDs(wellID));
                scheduleJS.control(k).W(wellID).val = 2;
            else
                % scheduleJS.control(k).W(wellID).val = -1 * well_fluid_rate_data(k, W_WellIDs(wellID));
                scheduleJS.control(k).W(wellID).val = -1e-2;
            end
            
            % scheduleJS.control(k).W.lims.bhp = 40e6;
            % scheduleJS.control(k).W.WI = scheduleJS.control(k).W.WI*factorWI;
        end
    end
    schedule = scheduleJS;
end

