function [state, bc, model, schedule, W, bestFaceIDs] = InitializeReservoir_SetBC_CreateWells_CreateSchedule(mainDir, model,  G, fluid,  rock)
    %{
        This function initialize the reservoir simulation for MRST GoM
        project. It also performs the following tasks:

        1) sets the boundary conditions
        2) Create the wells and their schedules
    %}

    %% Creating initial pressure and saving it
    g = norm(gravity);
    water_column = 200;
    p_r = g*fluid.rhoWS*water_column ;
    z_0 = 0; z_max = max(abs(G.cells.centroids(:,3)));
    equil  = ode23(@(z,p) g .* fluid.bW(p)*fluid.rhoWS, [z_0, z_max], p_r);
    p0 = reshape(deval(equil, abs(G.cells.centroids(:,3))), [], 1);  
    %%%p0=1000*g*G.cells.centroids(:,3);
    
    %% Plotting initial pressure 
    figure(ceil(rand)*100)
    clf()
    plot(p0/1e6, G.cells.centroids(:,3), 'o')
    hold on
    plot(fluid.rhoWS*g*abs(G.cells.centroids(:,3))/1e6 + p_r/1e6, (G.cells.centroids(:,3)),'r-', 'LineWidth',3)
    title('Initial pressure')
    xlabel('Initial pressure (MPa)')
    ylabel('Depth (m)')

    
    %% Saving initial reservoir pressures for later
    data.pressure=p0;
    data.porosity=G.cells.cellPoro;
    data.s=zeros(length(data.pressure),2);
    data.s(:,1)=1;
    data.t=0;
    data.FlowProps.ComponentPhaseDensity{1}=fluid.bW(p0)*fluid.rhoWS; % water density
    data.FlowProps.ComponentPhaseDensity{2}=[]; % water density  %% This is necessary otherwise the matlab julia read will throw an error
    data.FlowProps.ComponentPhaseDensity{3}=[]; % water density  %% This is necessary otherwise the matlab julia read will throw an error
    data.FlowProps.ComponentPhaseDensity{4}=fluid.bO(p0)*fluid.rhoOS; % CO2 density
    
    numChanged = 0;
    % modify state0 as discussed 
%     for curCell = 1 : size(G.cells.centroids(:, 3), 1)
%         if G.cells.centroids(curCell, 3) > -2000
%             data.s(curCell, 1) = 0;
%             data.s(curCell, 2) = 1;
%             numChanged = numChanged + 1;
%         end
%     end

    G_cells_picked = any(G.cells.centroids(:, 3) > -2000, 2);
    data_cells_picked = any(data.s(:, 2) == 1, 2);
    correct_cells = (G_cells_picked == data_cells_picked);
    
    dirToSave=strcat(mainDir, 'fromMatlab/mrstResults/data/');
    if (isfolder(dirToSave) == 0)
       mkdir(fullfile(mainDir, 'fromMatlab', 'mrstResults', 'data'));
    end
    fileToSave=strcat(dirToSave, 'state0')
    save(fileToSave,'data');

    
    
    %% Setting the boundary conditions
    indFaces = find(G.faces.neighbors(:,2) == 0);  %% Global index of the faces on the boundary of the domain
    indCell = G.faces.neighbors(indFaces,1);   %% cell index on the boundary of the domain

    %% now I set the boundary conditions as being pressure 
    bc=[];
    %s0 = zeros(length(indCell),2);
    %s0(:,1) = 1;
    %bc=addBC(bc, indFaces, 'pressure', p0(indCell), 'sat', s0 );
    %{
    figure(1)
    clf()
    scatter3(G.faces.centroids(ind, 1), G.faces.centroids(ind,2), G.faces.centroids(ind,3),'s')
    hold on
    scatter(G.cells.centroids(indCell, 1), G.cells.centroids(indCell, 2), 'o', '*k' )
    xlim([min(G.faces.centroids(:,1)), max(G.faces.centroids(:,1))] )
    ylim([min(G.faces.centroids(:,2)), max(G.faces.centroids(:,2))] )
    %}
    
    %% Create well and schedule here 
    % [schedule, W, tempList] = Create_Well_and_SimulationSchedule(mainDir, G, bc, rock);
  
    
    
    %% Initial reservoir state
    state = initResSol(G, p0, data.s );
    %state.flux = zeros(G.faces.num,2);
    state.wellSol = initWellSolAD(W, model, state);
    
    %% Setting the pore volume of the cells near the boundary to a very large value
    %model.operators.pv(indCell) = model.operators.pv(indCell)*1e4;
    

end