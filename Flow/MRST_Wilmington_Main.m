
%% MRST setup for the Exxon GoM small volume testing
%  1) multiphase flow
%  2) Immiscible flow

close all; clear all; clc

%%Gettting hostname
hostNameJS=char(java.net.InetAddress.getLocalHost.getHostName)
clusterHostName=split(hostNameJS,'.')

%% Supress all figures in the case of the cluster
%{
if strcmp(clusterHostName{end}, 'cluster')  || contains(hostNameJS, "node") %% slurm
    set(groot,'defaultFigureVisible','off')
end
%}

%% Computing initial pressure for the global domain
%ComputeInitialPressureGlobalDomain

%% Starting main function here
close all; clear all; clc 

%%Gettting hostname
hostNameJS=char(java.net.InetAddress.getLocalHost.getHostName)
clusterHostName=split(hostNameJS,'.')

%% Supress all figures in the case of the cluster
%{
if strcmp(clusterHostName{end}, 'cluster')  || contains(hostNameJS, "node") %% slurm
    set(groot,'defaultFigureVisible','off')
end
%}


% %% mainDir where everything will be saved at
mainDir=strcat(pwd(),'/')


%% Load MRST software libraries
run /Users/annikahuprikar/Downloads/mrst-2021a/startup.m

    
%% Directory containing induced seismicitiy libraries (Josimar)
%addpath('/home/josimar/Documents/Work/Projects/InducedSeismicity/LIB/InducedSeismicityTools/MATLAB/LIB/MATLAB');
%addpath('/home/josimar/Documents/Work/Projects/InducedSeismicity/LIB/InducedSeismicityTools/MATLAB/LIB/MATLABfunctions');

%%
mrstModule add ad-core ad-blackoil ad-props linearsolvers


%% Always turn on gravity first 
gravity reset on;


%% Directory to save MRST results, meaning the states (pressure, saturation etc)
dataFolderOutputSimulationState=strcat(mainDir, 'fromMatlab/mrstResults/');
if ~exist(dataFolderOutputSimulationState, 'dir')
    mkdir(dataFolderOutputSimulationState)
end

dirFigures=strcat(mainDir, 'Figures/');
if ~exist(dirFigures, 'dir')
    mkdir(dirFigures)
end

fprintf('\n \t Directory to output MRST results = %s \n', dataFolderOutputSimulationState)

%% Create handler object. Right now it does not output or store anything in memory because I do it myself in the loop
handler = ResultHandler('writeToDisk', true, 'storeInMemory', true, 'dataFolder', 'data', 'dataDirectory', dataFolderOutputSimulationState, 'cleardir', true);

%% Directory containing grid structure from Julia
JuliaFile=strcat(pwd,'/fromJulia/CellStruct.mat');
G=JuliaMRSTmapping(JuliaFile);
G = computeGeometry(G); % DONT USE WITH FAULTS
G.cells.centroids(:, 3) = abs(G.cells.centroids(:, 3));
G.faces.centroids(:, 3) = abs(G.faces.centroids(:, 3));
G.nodes.coords(:, 3) = abs(G.nodes.coords(:, 3));
fprintf('\n \t Julia file containing grid structure = %s \n ', JuliaFile)


%% Saving the cell centroid values for later QCing if needed
centroid=G.cells.centroids;
dirToSave=strcat(mainDir, 'fromMatlab/mrstResults/data/');
if (isfolder(dirToSave) == 0)
   mkdir(fullfile(mainDir, 'fromMatlab', 'mrstResults', 'data'));
else
    dirTmp=strcat(mainDir, 'fromMatlab/mrstResults/data/');
   % tmp=join(['rm ' dirTmp '*mat']);
   % system(tmp);
end
fileToSave=strcat(dirToSave, 'cellCentroid')
save(fileToSave,'centroid');


%% Reading fluid properties and creating fluid model object here
%inputFile=strcat(mainDir,'inputData/fluidProperties/GoM_FluiProperties.DATA')
%fluidOLD = ReadFluidProperties_GoM(mainDir, inputFile);

%inputFile=strcat(mainDir,'inputData/fluidProperties/ftest_co2brine.DATA')
%inputFile=strcat(mainDir,'inputData/fluidProperties/josimarGoM.DATA.txt')

% inputFile=strcat(mainDir,'inputData/fluidProperties/GoM_NoDissolution_With_Pc.txt') % INSERT NEW FILE 3/22/22
inputFile=strcat(mainDir,'Wilmington_MultiphaseFlow.dat');
fluid = ReadFluidProperties_Wilmington(mainDir, inputFile);

%%
%{
sw=linspace(0,1,100);
figure(100)
clf()
subplot(121)
plot(sw, fluidOLD.krW(sw),'-r')
hold on
plot(sw, fluid.krW{1}(sw),'--b')

subplot(122)
plot(sw, fluidOLD.krOW(1 - sw),'-r')
hold on
plot(sw, fluid.krOW{1}(1 - sw),'--b')
%}

%%

%% Define rock model
% rock = makeRock(G, G.cells.cellPerm(:,1), G.cells.cellPoro);
rock = makeRock(G, [20 20 4].*milli*darcy, 0.15); % CHANGE ROCK PROPERTIES LATER 3/22/22
% permeability = 200 mD, [200, 200, 40]
% porosity = 0.15
%rock = makeRock(G, G.cells.cellPerm(:,1), 0.1);
%rock = makeRock(G,1e-13, 0.1);

%% Fixing the capillary pressure here
%% Don't use this for the ISGS project because the matrix perm are very low.
%[fluid, rock] = FixFaultCapillaryPressure(G, rock, fluid);


%% Create model here Set up model and initial state.
model = TwoPhaseOilWaterModel(G, rock, fluid); % TWO PHASE FLOW
% model = GenericBlackOilModel(G, rock, fluid, 'gas', false, 'gravity', [0 0 -norm(gravity)]); % SINGLE PHASE FLOW
% model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false);
% model.oil=0;

%% Set up mex-accelerated backend and reduced variable set for wells
% We use the Diagonal autodiff backend to calculate derivatives. By
% default, this uses only Matlab, so we also set the "useMex" flag to be
% true. It will then use C++ versions of most discrete operators during
% assembly.
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
model = model.validateModel();
%model.FacilityModel.primaryVariableSet = 'val';
model.OutputStateFunctions = {'Mobility', 'ComponentTotalMass','RelativePermeability', 'ComponentPhaseDensity'};
ncomp = 2;
%ncomp = model.getNumberOfComponents();
solver = AMGCL_CPRSolverAD('tolerance', 1e-4, 'block_size', ncomp, ...
                            'useSYMRCMOrdering', true, ...
                            'coarsening', 'aggregation', 'relaxation', 'ilu0', ...
                            'maxIterations', 200);

nls = NonLinearSolver('LinearSolver', solver);
%model.toleranceMB = 1e-4;
%model.verbose = true;
%% Run assembly tests (will compile backends if required)
% Note: Benchmarks are obviously not accurate if compilation is performed.
%testMexDiagonalOperators(model, 'block_size', 3);

%% Initialize reservoir state and set boundary conditions, create wells and schedule
% [state0, bc, model, schedule, W, bestFaceIDs] = InitializeReservoir_SetBC_CreateWells_CreateSchedule(mainDir, model, G, fluid, rock);
bc = []
state0 = FunctionComputeInitialPressure_MultiphaseFlow(mainDir, model, G, fluid, 2200, 5, 5);
[schedule, W, bestFaceIDs] = addWellsAndSimulationSchedule(mainDir, G, bc, rock);
model
totalNumberOfCells = length(G.cells.cellPerm)
%% Solve problem
%nls = NonLinearSolver('useLinesearch', true);
%[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'nonlinearsolver', nls, 'outputHandler', handler);
%[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'outputHandler', handler);

%% Set up a compiled linear solver
% We use a AMGCL-based CPR solver to solve the problem. It is sufficient to
% have a working C++ compiler, the AMGCL repository and BOOST available to
% compile it. See the documentation of amgcl_matlab for more details on how
% to set up these paths.


%% Run a parallel simulation with four threads
maxNumCompThreads(5);
%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls, 'outputHandler', handler);

% [~, ~, report] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls, 'outputHandler', handler);
[wellSols, states, report] = simulateScheduleAD(state0, model, scheduleJS, 'outputHandler', handler, 'OutputMinisteps', true, 'NonLinearSolver', nls);


% STEPS/IDEAS:
% states: saturation and pressure in every cell
% report: model simulations

%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, 'outputHandler', handler);

% wellSols.qOs - analogous, for oil
% wellSols.qWs - one val / timestep of how much well produced of water

% see this overtime

% sigma production of oil and water for EVERY WELL
% total oil production of wells vs. total oil production of field 
% analogous for water



% make plots starting with one well (producer)
% play with oil-water contact - how do results change?

% repeat with more producers (larger magnitude)


% set timespan to 30 years (when only production was taking place)

% sanity check: run one injector to make sure behavior is as expected

%% Export results for later comparison with Julia
simTime=report.ReservoirTime;
FileToSaveData=strcat(dataFolderOutputSimulationState,'data/simTime.mat');
save(FileToSaveData, 'simTime')


%% Show Overall Difference in Pressure from First to Last Time Step
load('/Users/annikahuprikar/Desktop/Spring2021/HSURV PRISE Summer 2021/Flow/fromMatlab/mrstResults/data/state0.mat')
data1=data;

load('/Users/annikahuprikar/Desktop/Spring2021/HSURV PRISE Summer 2021/Flow/fromMatlab/mrstResults/data/state20.mat')
data2=data;

figure()
clf()
plotCellData(G, (data2.pressure - data1.pressure)/1e6)

colorbar
view(3)
ax = gca;
ax.ZDir = 'normal';

hold on
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
hold off

diffFig = figure;
scatter((data2.pressure - data1.pressure) / 1e6, G.cells.centroids(:, 3))
title("Change in Pressure")
xlabel("Pressure Difference (MPa)")
ylabel("Depth")

cell_p_sat = analyzeCell(mainDir, 3658);
cell_neighbors = getCellNeighbors(G, 3658);
figure()
scatter(cell_p_sat(:, 1), cell_p_sat(:, 2))
scatter(cell_p_sat(:, 1), cell_p_sat(:, 3))

for entry = 1 : numel(cell_neighbors) 
    cell_p_sat = analyzeCell(mainDir, cell_neighbors(entry));
    figure()
    scatter(cell_p_sat(:, 1), cell_p_sat(:, 2))
    scatter(cell_p_sat(:, 1), cell_p_sat(:, 3))
end