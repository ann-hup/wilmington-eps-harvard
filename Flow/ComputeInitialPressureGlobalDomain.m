%% Computes the initial pore pressure for the PyLith one way coupled domain using the formation volume factors 

% To ensure consitency between the global and the filtered domains, I have
% to make sure that the initial pressures for the PyLith one way coupled
% model in the global domain (not only in the filtered domain) are computed
% based on the formation volume factors

%%
close all; clear all; clc 

%%Gettting hostname
hostNameJS=char(java.net.InetAddress.getLocalHost.getHostName)
clusterHostName=split(hostNameJS,'.')

%% mainDir where everything will be saved at
mainDir=strcat(pwd(),'/')

%% Directory to save MRST results, meaning the states (pressure, saturation etc)
dataFolderOutputSimulationState=strcat(mainDir, 'fromMatlab/mrstResults/');

if ~exist(dataFolderOutputSimulationState, 'dir')
    mkdir(dataFolderOutputSimulationState)
end

dirFigures=strcat(mainDir, 'Figures/');

if ~exist(dirFigures, 'dir')
    mkdir(dirFigures)
end


if strcmp(hostNameJS , 'Josimars-iMac.local')
     
    %% Start MRST
    run /Users/josimar/Documents/Work/Projects/Induced_Seismicity/Development/MRST_Poroelasticity/LIB/MRST_Dev/mrst-2020a/startup.m
    
    %% Directory containing induced seismicitiy libraries (Josimar)
    addpath('/Users/josimar/Documents/Work/Projects/Induced_Seismicity/LIB/InducedSeismicityTools/MATLAB/LIB/MATLAB')
    addpath('/Users/josimar/Documents/Work/Projects/Induced_Seismicity/LIB/InducedSeismicityTools/MATLAB/LIB/MATLABfunctions');

elseif strcmp(clusterHostName{end}, 'cluster')  || contains(hostNameJS, "node") %% slurm
   
    %% Start MRST
    run /nobackup1/josimar/Software/ReservoirSimulation/MRST/mrst-2020a/startup.m
    
    %% Directory containing induced seismicitiy libraries (Josimar)
    addpath('/nobackup1/josimar/Projects/InducedSeismicity/LIB/InducedSeismicityTools/MATLAB/LIB/MATLAB')
    addpath('/nobackup1/josimar/Projects/InducedSeismicity/LIB/InducedSeismicityTools/MATLAB/LIB/MATLABfunctions');

elseif strcmp(hostNameJS , 'GRS-SERVER')
    
    run C:\Users\josimar\Documents\Work\Projects\InducedSeismicity\Resources\Software\MRST_GRS_Server\mrst-2020a\startup.m
    
    %% Directory containing induced seismicitiy libraries (Josimar)
    addpath('C:\Users\josimar\Documents\Work\Projects\InducedSeismicity\Resources\LIB\InducedSeismicityTools\MATLAB\LIB\MATLAB')
    
    %% Directory containing libraries related to the GoM work. Some of them can be useful for other projects (e.g ISGS)
    addpath('C:\Users\josimar\Documents\Work\Projects\InducedSeismicity\Resources\LIB\InducedSeismicityTools\MATLAB\LIB\MATLABfunctions')

elseif strcmp(hostNameJS , 'epsgrs5') %% epsgrs5
    %% Start MRST
    run /home/josimar/Documents/Work/Resources/Software/ReservoirSimulation/MRST/MRST_OWC_only/mrst-2020b/startup.m
    
    %% Directory containing induced seismicitiy libraries (Josimar)
    addpath('/home/josimar/Documents/Work/Projects/InducedSeismicity/LIB/InducedSeismicityTools/MATLAB/LIB/MATLAB');
    addpath('/home/josimar/Documents/Work/Projects/InducedSeismicity/LIB/InducedSeismicityTools/MATLAB/LIB/MATLABfunctions');
    
end %%if

%% Start MRST

mrstModule add ad-core ad-blackoil ad-props linearsolvers

%% Always turn on gravity first 
gravity reset on;


%% Reading fluid properties and creating fluid model object here
%inputFile=strcat(mainDir,'inputData/fluidProperties/GoM_FluiProperties.DATA')
%inputFile=strcat(mainDir,'inputData/fluidProperties/BENCH_SPE1.DATA')
%inputFile=strcat(mainDir,'inputData/fluidProperties/ftest_co2brine.DATA')
%fluid = ReadFluidProperties_GoM(mainDir, inputFile);

inputFile=strcat(mainDir,'inputData/fluidProperties/GoM_NoDissolution_With_Pc.txt')
fluid = ReadFluidProperties_GoM_withoutDissolution_WithPc(mainDir, inputFile);

%% Read cell centroids for the global domain
fileName=strcat(mainDir,'fromJulia/cellCentroidsGlobalDomain.mat');
load(fileName);
cellCentroids = abs(cellCentroids);

%% Creating initial pressure and saving it
g = norm(gravity);
water_column = 200;
p_r = g*fluid.rhoWS*water_column ;
z_0 = 0; z_max = max(cellCentroids (:,4));
equil  = ode23(@(z,p) g .* fluid.bW(p)*fluid.rhoWS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, cellCentroids (:,4)), [], 1);  


fileNameToSavePressure=strcat(mainDir,'fromMatlab/mrstResults/InitialPressureGlobalDomain.mat');
save(fileNameToSavePressure,'p0', 'cellCentroids')

%set(0, 'DefaultLineLineWidth', 2);
%set(0,'defaultAxesFontSize',16)
%set(groot,'defaultFigurePaperPositionMode','manual')

figureName=strcat(mainDir,'Figures/InitialPressureGlobalDomain.png');

stepPlot=1:100:size(cellCentroids,1);

figure(10)
clf()
scatter(p0(stepPlot)/1e6, cellCentroids(stepPlot,4), 'bo')
hold on
plot(1070*9.8*cellCentroids(stepPlot,4)/1e6, cellCentroids(stepPlot,4), '-r')
title('Initial pressure for the global domain')
xlabel('P0 (MPa)')
ylabel('Depth (m)')
legend('Using PVT', 'Using rho=1000 kg/m3', 'Location', 'northwest')

saveas(gcf,figureName)



%% Now computes the initial fluid density in the global domain
rhoW0 = fluid.bW(p0)*fluid.rhoWS;

fileNameToSavePressure=strcat(mainDir,'fromMatlab/mrstResults/InitialFluidDensityGlobalDomain.mat');
save(fileNameToSavePressure,'rhoW0', 'cellCentroids')

figureName=strcat(mainDir,'Figures/InitialFluidDensityGlobalDomain.png');
figure(11)
clf()
scatter(rhoW0(stepPlot), cellCentroids(stepPlot,4), 'bo')
hold on
%plot(1070*9.8*cellCentroids(stepPlot,4)/1e6, cellCentroids(stepPlot,4), '-r')
title('Initial fluid density for the global domain')
xlabel('water density (g/cm3)')
ylabel('Depth (m)')
%legend('Using PVT', 'Using rho=1000 kg/m3', 'Location', 'northwest')

saveas(gcf,figureName)


