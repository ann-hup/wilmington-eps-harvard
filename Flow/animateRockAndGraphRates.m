function outputFig = animateRockAndGraphRates(mainDir, G)
%ANIMATEROCK Summary of this function goes here
%   Detailed explanation goes here

    well_cell_data = dlmread(strcat(mainDir, 'superwellData/well_cell_mapping.dat'));  
    well_fluid_rate_data = dlmread(strcat(mainDir, 'superwellData/superwell_fluid_rates.dat')); 
    timestep_info = dlmread(strcat(mainDir, 'dataExported/Injection/date.dat')); 
    
    injector_locs = any(well_cell_data == 1, 2);
    producer_locs = ~injector_locs;
    
    injector_indexes = well_cell_data(injector_locs, 1);
    producer_indexes = well_cell_data(producer_locs, 1);
    
    injector_rates_only = well_fluid_rate_data(:, injector_indexes);
    producer_rates_only = well_fluid_rate_data(:, producer_indexes);
    
    injector_sum_dt = sum(injector_rates_only, 2);
    producer_sum_dt = sum(producer_rates_only, 2);
    diff_sum_dt = injector_sum_dt - producer_sum_dt;
    
    loc_1936 = find(timestep_info == 1936); loc_1975 = find(timestep_info == 1975); 
    loc_1985 = find(timestep_info == 1985);
    
    statesDir = strcat(mainDir, 'fromMatlab/mrstResults/data/');
    datesFile = dlmread(strcat(mainDir, 'dataExported/Injection/date.dat'));
    initialMatObj = matfile(strcat(statesDir, 'state0.mat'));
    initialMatObjStruct = initialMatObj.data(1, :);
    az = 5; el = 20; increase_el = 1; txt = '';
    
    rgb = [ ...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   191
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;

    rgb2 = [ ...
    110 38 14
    136    8   8
    196   30   58
   210   43   43
   236   88   0
   227 115 94
   250   128   114
   250   213   165
   255   248   220
   255 255 255 % 10
   236 255 220
   193 225 193
   147 197 114
   64 181 173
   0 0 255
   63 0 255
   0 71 171
   127 0 255
   93 63 211
   0 0 0 ] / 255;
    
    outputFig = figure('Name', 'Rock Animation');
    subplot(2, 1, 1);
    plotCellData(G, (initialMatObjStruct.pressure)/1e6)
    hold on
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    ax = gca;
    ax.ZDir = 'normal';
    % caxis([-20 10])
    % colorbar('Limits', [-20 8])
    % OR 
    %{
    colormap(rgb)
    %h = colorbar
    %h.Ticks = linspace(-10, 10, 20)
    colorbar
    set(gca,'CLim',[-5 20])
    %}
    colorbar
    %set(gca,'CLim',[-5 25])
    view(az, el)
    
    subplot(2, 1, 2)
    plot(timestep_info, injector_sum_dt);
    hold on
    plot(timestep_info, producer_sum_dt);
    title('Total Injection and Production Rates in m^3/s');
    legend('injection', 'production', 'Location', 'northwest');
    persistent lineHandle 
    if isempty(lineHandle) || ~isvalid(lineHandle)
        lineHandle = xline(1936, 'r-'); 
    end
    hold off
    
    
    outputFig.Position(1:4) = [300, 300, 1400, 1100];

    for stateK = 1:1008
        subplot(2, 1, 1);
        stateNum = ['state', int2str(stateK), '.mat'];
        matObj = matfile(strcat(statesDir, stateNum));
        matObjStruct = matObj.data(1, :);

        plotCellData(G, (matObjStruct.pressure - initialMatObjStruct.pressure)/1e6)
        hold on
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        ax = gca;
        ax.ZDir = 'normal';
        %caxis([-20 10])
        %colorbar('Limits', [-20 8])
        % OR 
        %{
        colormap(rgb)
        %h = colorbar
        %h.Ticks = linspace(-10, 10, 20)
        colorbar
        set(gca,'CLim',[-5 20])
        %}
        colorbar
        %set(gca,'CLim',[-5 25])
        
        if mod(1 + stateK, 12) == 1
            txt = datesFile(stateK);
            disp(txt); 
        end
        
        text(4.02e5, 3.73e6, 300, int2str(txt), 'color','black','fontsize', 24);
        
        if stateK > 560 && stateK < 710
            az = mod(az + 2, 360);
        end
        
        if stateK > 350 && stateK < 430
           if el >= 80
               increase_el = 0;
           end
           if el <= 0
               increase_el = 1;
           end
           
           if increase_el == 1
               el = el + 2;
           else
               el = el - 2;
           end
        end
        
        view(az, el)
        outputFig.Position(1:4) = [300, 300, 1400, 1100];

        pause(0.01)
        disp(stateK)
 
        if stateK ~= 1008
            hold on
            cla
        end
        
%         if datesFile(stateK) == 1936
%             subplot(2, 1, 2);
%             plot(timestep_info(loc_1936 : loc_1975), injector_sum_dt(loc_1936 : loc_1975));
%             hold on
%             plot(timestep_info(loc_1936 : loc_1975), producer_sum_dt(loc_1936 : loc_1975));
%             title('1936-1975 Total Injection and Production Rates in m^3/s');
%             legend('injection', 'production', 'Location', 'northwest');
%             if isempty(lineHandle) || ~isvalid(lineHandle)
%                 lineHandle = xline(1936 + (1/12), 'r-'); 
%             end
%             hold off
%         end
%         if datesFile(stateK) == 1975
%             subplot(2, 1, 2);
%             plot(timestep_info(loc_1975 : loc_1985), injector_sum_dt(loc_1975 : loc_1985));
%             hold on
%             plot(timestep_info(loc_1975 : loc_1985), producer_sum_dt(loc_1975 : loc_1985));
%             title('1975-1985 Total Injection and Production Rates in m^3/s');
%             legend('injection', 'production', 'Location', 'northwest');
%             if isempty(lineHandle) || ~isvalid(lineHandle)
%                 lineHandle = xline(1975 + (1/12), 'r-'); 
%             end
%             hold off
%         end
%         if datesFile(stateK) == 1985
%             subplot(2, 1, 2);
%             plot(timestep_info(loc_1985 : end), injector_sum_dt(loc_1985 : end));
%             hold on
%             plot(timestep_info(loc_1985 : end), producer_sum_dt(loc_1985 : end));
%             title('1985-2020 Total Injection and Production Rates in m^3/s');
%             legend('injection', 'production', 'Location', 'northwest');
%             if isempty(lineHandle) || ~isvalid(lineHandle)
%                 lineHandle = xline(1985 + (1/12), 'r-'); 
%             end
%             hold off
%         end
%         if (datesFile(stateK) > 1936 && datesFile(stateK) < 1975) || ...
%                 (datesFile(stateK) > 1975 && datesFile(stateK) < 1985) || ...
%                 (datesFile(stateK) > 1985)
%             subplot(2, 1, 2);
%             lineHandle.Value = lineHandle.Value + (1/12);
%         end
    end
end


%{
REMINDER
- have mrstResults zip file
- have data in well_cell_data and well_fluids_rates files done
- shouldn't need to run this simulation again ...

TO DO
- perhaps bestFaces visual
- pick 5 cells to model pressure changes in (internal cells)
- work on presentation

FALL 2021: Multiphase Flow (oil and water)
- oil fluid properties and reservoir properties - better understand
- viscosity, formation volume factor, density of oil - how param change
with pressure (value API 26 - viscosity)
- literature
- read chapter in textbook
%}