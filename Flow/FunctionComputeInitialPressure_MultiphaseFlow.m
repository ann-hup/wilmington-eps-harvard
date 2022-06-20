function state0 = FunctionComputeInitialPressure_MultiphaseFlow(mainDir, model, G, fluid, depthOWC, matIdReservoir, matIdToPlot)

    dirToSave=strcat(mainDir,'fromMatlab/mrstResults');
    if (isfolder(dirToSave) == 0)
       mkdir(fullfile(mainDir, 'fromMatlab', 'mrstResults'));
    end
    
    if (isfolder(strcat(mainDir,'Figures')) == 0)
       mkdir(fullfile(mainDir, 'fromMatlab', 'Figures'));
    end
    
    if depthOWC < 0
        error("Oil water contact must be negative!!")
    end
    
    %% Getting cells Ids for the reservoir
    idToPlotResults=[];
    for i=1:length(matIdToPlot)
        idTmp = find(G.cells.materialId == matIdToPlot(i));    
        if isempty(idTmp) == 0
            idToPlotResults=[idToPlotResults; idTmp];
        end
    end
    
    %% Initialize model for the multiphase flow case
    %idReservoir = find(G.cells.materialId <= 4);
    %idNonReservoir = find(G.cells.materialId > 4);
    
    
    matIdNonReservoir=unique(G.cells.materialId);
    for i=1:length(matIdReservoir)
        matIdNonReservoir(matIdNonReservoir == matIdReservoir(i)) = [];
    end

   
    %% Getting cells Ids for the reservoir
    idReservoir=[];
    for i=1:length(matIdReservoir)
        idTmp = find(G.cells.materialId == matIdReservoir(i));    
        if isempty(idTmp) == 0
            idReservoir=[idReservoir; idTmp];
        end
    end
    
    idReservoirOriginal=idReservoir;

    %%
    
    %% Getting cells Ids for the non-reservoir part
    idNonReservoir=[];
    for i=1:length(matIdNonReservoir)
        idTmp = find(G.cells.materialId == matIdNonReservoir(i));    
        if isempty(idTmp) == 0
            idNonReservoir=[idNonReservoir; idTmp];
        end
    end
    
    
    
    %% The method below is not very good
    %{
    %% Getting cells Ids for the reservoir
    idReservoir=[];
    for i=1:length(matIdReservoir)
        idTmp = find(G.cells.materialId == matIdReservoir(i));    
        if isempty(idTmp) == 0
            idReservoir=[idReservoir; idTmp];
        end
    end
    
    idReservoirOriginal = idReservoir;
    
    %% Now I have to remove the basement cells that are above the OWC contact
    idToKeep=zeros(length(idReservoir),1);
    countIdToKeep=1;
    for i=1:length(idReservoir)
        if (G.cells.materialId(i) == matIdBasement) %% basement below ranger
            if G.cells.centroids(i,3) > depthOWC
                %idToKeep = [idToKeep; i];
                idToKeep(countIdToKeep) = i;
                countIdToKeep = countIdToKeep + 1;
            end
        else
            %idToKeep = [idToKeep; i];
            idToKeep(countIdToKeep) = i;
            countIdToKeep = countIdToKeep + 1;

        end
    end

    idToKeep = idToKeep(1:countIdToKeep - 1);
    idReservoir = idReservoir(idToKeep);

    idAll=1:G.cells.num;
    idNonReservoir = setdiff(idAll, idReservoir);
    %}
    
    %%

    tmp = 1:model.G.cells.num;
    tmp= tmp(idReservoir);
    ddepth= 0;
    p_ref = 0;
    region = getInitializationRegionsBlackOil(model, depthOWC, 'datum_depth', ddepth, 'datum_pressure', p_ref, 'cells', tmp );
    state1 = initStateBlackOilAD(model, region);

%     tmp = 1:model.G.cells.num;
%     tmp= tmp(idNonReservoir);
%     ddepth=0;
%     p_ref = 0;
%     region = getInitializationRegionsBlackOil(model, depthOWC*0, 'datum_depth', ddepth, 'datum_pressure', p_ref, 'cells', tmp );
%     state2 = initStateBlackOilAD(model, region);

    %state0.pressure = state1.pressure + state2.pressure;
    %state0.s = state1.s + state2.s;

    state0.pressure = zeros(G.cells.num,1);
    state0.s = ones(G.cells.num,2);

    state0.pressure(idReservoir) = state1.pressure(idReservoir);
%     state0.pressure(idNonReservoir) = state2.pressure(idNonReservoir);
    state0.s(idReservoir,:) = state1.s(idReservoir,:);
%     state0.s(idNonReservoir,:) = state2.s(idNonReservoir,:);
    
    
 

    %{
    idReservoir = find(G.cells.materialId < 3);
    tmp = 1:model.G.cells.num;
    ddepth= 0;
    p_ref = 0;
    region = getInitializationRegionsBlackOil(model, depthOWC, 'datum_depth', ddepth, 'datum_pressure', p_ref, 'cells', tmp );
    state1 = initStateBlackOilAD(model, region);
    state0 = state1;
    %}

    %% Saving results now
    cellCentroids = G.cells.centroids;
    s0 = state0.s;
    p0 = state0.pressure;
    rhoW0 = fluid.bW(p0) .* fluid.rhoWS .* state0.s(:,1) + fluid.bO(p0) .* fluid.rhoOS .* state0.s(:,2)  ; %% Total fluid density
    
    fileNameToSavePressure=strcat(mainDir,'fromMatlab/mrstResults/InitialPressureGlobalDomain.mat');
    save(fileNameToSavePressure,'p0', 'cellCentroids')

    fileNameToSavePressure=strcat(mainDir,'fromMatlab/mrstResults/InitialFluidDensityGlobalDomain.mat');
    save(fileNameToSavePressure,'rhoW0', 'cellCentroids')
    
    fileNameToSavePressure=strcat(mainDir,'fromMatlab/mrstResults/InitialSaturationGlobalDomain.mat');
    save(fileNameToSavePressure,'s0', 'cellCentroids')

    
    %% Saving data for visualization
    data.pressure=p0;
    data.porosity=G.cells.cellPoro;
    data.s=state0.s;
    data.FlowProps.Density = rhoW0  ; %% Total fluid density
    data.t=0;

    dirToSave=strcat(mainDir, 'fromMatlab/mrstResults/data/');
    if (isfolder(dirToSave) == 0)
       mkdir(fullfile(mainDir, 'fromMatlab', 'mrstResults', 'data'));
    end
    fileToSave=strcat(dirToSave, 'state0')
    save(fileToSave,'data');

    
    
    %% Plotting 
    %set(0, 'DefaultLineLineWidth', 2);
    %set(0,'defaultAxesFontSize',16)
    %set(groot,'defaultFigurePaperPositionMode','manual')

    
    figureName=strcat(mainDir,'Figures/InitialPressureGlobalDomain.png');

    stepPlot=1:100:G.cells.num;

    figure(10)
    clf()
    scatter(p0(stepPlot)/1e6, abs(G.cells.centroids(stepPlot,3)), 'bo')
    hold on
    plot(fluid.rhoWS*9.8*abs(cellCentroids(stepPlot,3))/1e6, abs(G.cells.centroids(stepPlot,3)), '-r')
    title('Initial pressure for the global domain')
    set(gca, 'YDir','reverse')
    xlabel('P0 (MPa)')
    ylabel('Depth (m)')
    legend('Using PVT', 'Using rho=1000 kg/m3', 'Location', 'east')

    saveas(gcf,figureName)


    figureName=strcat(mainDir,'Figures/InitialFluidDensityGlobalDomain.png');
    figure(11)
    clf()
    scatter(rhoW0(stepPlot), abs(G.cells.centroids(stepPlot,3)), 'bo')
    %hold on
    %plot(1070*9.8*cellCentroids(stepPlot,4)/1e6, cellCentroids(stepPlot,4), '-r')
    set(gca, 'YDir','reverse')
    title('Initial fluid density for the global domain')
    xlabel('water density (g/cm3)')
    ylabel('Depth (m)')
    %legend('Using PVT', 'Using rho=1000 kg/m3', 'Location', 'northwest')

    saveas(gcf,figureName)
    
    G.nodes.coords(:,3) = abs(G.nodes.coords(:,3));
    figureName=strcat(mainDir,'Figures/OWC.png');
    figure(12)
    clf()
    plotCellData(G, data.s(:,2), idReservoirOriginal)
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    %set(gca, 'ZDir','reverse')
    title('Initial oil saturation (Reservoir Only)')
    colorbar
    view(3)
    saveas(gcf,figureName)


    figureName=strcat(mainDir,'Figures/OWC_entireDomain.png');
    figure(13)
    clf()
    plotCellData(G, data.s(:,2))
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    %set(gca, 'ZDir','reverse')
    title('Initial oil saturation (Entire domain)')
    colorbar
    view(3)
    saveas(gcf,figureName)
    
    figureName=strcat(mainDir,'Figures/OWC_SelectedZones.png');
    figure(14)
    clf()
    plotCellData(G, data.s(:,2), idToPlotResults)
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    %set(gca, 'ZDir','reverse')
    title('Initial oil saturation (Selected zones)')
    colorbar
    view(3)
    saveas(gcf,figureName)

    
    
end
