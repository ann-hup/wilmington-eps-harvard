function fluid = ReadFluidProperties_Wilmington(mainDir, inputFile, plotFigures)
    %{
        Read and plot the fluid properties for the GoM simulations
    %}

    
    set(0, 'DefaultLineLineWidth', 2);
    set(0,'defaultAxesFontSize',16)
    set(groot,'defaultFigurePaperPositionMode','manual')

    %%
    mrstModule add ad-props deckformat mrst-gui ad-core ad-blackoil
    mrstModule add deckformat ad-props

    fileName=strcat(inputFile)
    deck = readEclipseDeck(fileName);
    deck = convertDeckUnits(deck);

    % Create a special ADI fluid which can produce differentiated fluid
    % properties.
    fluid = initDeckADIFluid(deck)
    
    %% Checking if the user wants to plot the fluid information or not
    if ~exist('plotFigures','var')
     % third parameter does not exist, so default it to something
      plotFigures = 0;
    end

    if plotFigures == 1
        
        P=linspace(0e6,170e6, 100);
        %sW=linspace(0,1, 100);
        sW=linspace(min(fluid.krPts.w(:)),max(fluid.krPts.w(:)),100);

        dirToSave=strcat(mainDir, 'Figures/');
        if (isfolder(dirToSave) == 0)
            mkdir(fullfile(mainDir, 'Figures'));
        else

        figureName=strcat(mainDir,'Figures/FluidProperties.pdf');

        fig=figure('Renderer', 'painters', 'Position', [100 100 1300 900]);
        clf()
        orient(fig,'landscape')
        subplot(221)
        plot(P/1e6, fluid.bW(P))
        %plot(P/1e6, fluid.bO(P,0,1))
        title('1/Bw')
        xlabel('Pressure (MPa)')
        ylabel('1/Bw (sm3/rm3)')

        subplot(222)
        semilogy(P/1e6, fluid.muW(P))
        %plot(P/1e6, fluid.muO(P))
        %plot(P/1e6, fluid.muO(P,0,1))
        title('Water viscosity ')
        ylabel('viscosity (Pa * s)')
        xlabel('Pressure (MPa)')

        subplot(223)
        %plot(P/1e6, fluid.bG(P))
        plot(P/1e6, fluid.bO(P))
        title('1/Bg')
        xlabel('Pressure (MPa)')
        ylabel('1/Bg (sm3/rm3)')

        subplot(224)
        %plot(P/1e6, fluid.muG(P))
        semilogy(P/1e6, fluid.muO(P))
        title('CO2 viscosity ')
        ylabel('viscosity (Pa * s)')
        xlabel('Pressure (MPa)')
        saveas(gcf,figureName)

       

        figureName=strcat(mainDir,'Figures/RelativePermeability_Manuscript.pdf');
        figure(21)
        clf()
        plot(sW,fluid.krW(sW), 'DisplayName','Brine', 'LineWidth',2)
        hold on
        plot(sW,fluid.krOW(1-sW), 'r', 'DisplayName','Oil', 'LineWidth',2)
        xlabel('Water saturation')
        ylabel('Relative permeability')
        legend
        xlim([0,1])
        ylim([0,1])
        saveas(gcf,figureName)   

        x1=1450*9.8*1000/1e6;
        x2=1550*9.8*1000/1e6 + 10;

        figureName=strcat(mainDir,'Figures/FluidDensityVariationWithPressure.pdf');
        figure()
        clf()
        plot(P/1e6, fluid.bW(P) * fluid.rhoWS, '-b')
        hold on
        plot(P/1e6, fluid.bO(P) * fluid.rhoOS, '-r')
        hold on
        %plot(ptmp/1e6, rhoOfit, '--k')
        hold on
        %fill([x1;x1; x2;x2], [-200; 1400; 1400; -200], 'black')
        %alpha(0.1)
        %xlim([0,40])
        ylim([0,1200])
        xlabel('Pressure (MPa)')
        ylabel('Density (kg/m3)')
        %legend('Brine', 'CO2 original', 'CO2 modified', 'Location','SouthEast')
        legend('Brine', 'CO2', 'Location','SouthEast')
        saveas(gcf,figureName)


        figureName=strcat(mainDir,'Figures/FluidViscosityVariationWithPressure.pdf');
        figure()
        clf()
        semilogy(P/1e6, fluid.muW(P), '-b')
        hold on
        semilogy(P/1e6, fluid.muO(P), '-r')
        %fill([x1;x1; x2;x2], [0; 1.2e-3; 1.2e-3; 0], 'black')
        %alpha(0.1)
        %plot(ptmp/1e6, rhoOfit, '--k')
        xlabel('Pressure (MPa)')
        ylabel('Viscosity (Pa * s)')
        legend('Brine', 'CO2', 'Location','Northwest')
        saveas(gcf,figureName)
    
    end
    
end