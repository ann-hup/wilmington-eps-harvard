function [fluid, rock] = FixFaultCapillaryPressure(G, rock, fluid)

    %% Including now the fault capillary pressure 
    idFaults = G.cells.materialId < 0; % i.e. top "seal" layer (= caprock)
    Faults_satRegNum = 2;

    rock.regions.saturation = ones(G.cells.num,1); 
    rock.regions.saturation(idFaults) = Faults_satRegNum;

    %%% Function from LLuis Salo
    refPc.val  = fluid.pcOW{Faults_satRegNum};       % you would use the id for fault cells
    refPc.poro = max(rock.poro(idFaults));  
    refPc.perm = max(rock.perm(idFaults));  

    fault.poro = rock.poro(idFaults);                % must be passed in ascending order of global fault cell id
    fault.perm = rock.perm(idFaults);                % "

    fault.satRegNum = 2;                       

    fluid = assignFaultPcSimple(fluid, fault, refPc);   

    sW=linspace(0,1,100);
    figure(13)
    clf
    plot(sW, fluid.pcOW{2}(sW)/1e6)
    xlabel('Wate saturation')
    ylabel('Fault Capillary pressure (MPa)')


    sW=linspace(0,1,100);
    figure(14)
    clf
    plot(sW, fluid.pcOW{1}(sW)/1e6)
    xlabel('Wate saturation')
    ylabel('Matrix Capillary pressure (MPa)')

    %{
    pAll = linspace(0,170e6,100);
    figure(15)
    clf
    plot(pAll/1e6, fluid.rsSat(pAll))
    xlabel('Pressure (MPa)')
    ylabel('Solution gas-oil ratio')
    %}


end