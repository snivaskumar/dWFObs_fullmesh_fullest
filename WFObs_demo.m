clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%             WIND FARM OBSERVER (WFObs) by B.M. Doekemeijer
%                 Delft University of Technology, 2017
%              Repo: https://github.com/Bartdoekemeijer/WFObs
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%
%%   Quick use:
%     1. Specify your script preferences in lines 58 - 73
%     3. Specify your simulation settings file in line  76
%      a. if you want (OPTIONAL!), you can create a simulation case 
%         manually, then have a look at the ./configurations folder
%          i.  Use  '/Setup_sensors/WFObs_s_sensors.m' to specify measurement locations
%          ii. All WFSim model/site information is contained in '../WFSim/bin/core/meshing.m'
%     3. Press start.
%
%%   Relevant input/output variables
%     - configName: name of the simulation case that is to be simulated.
%     All simulation scenarios can be found in the '/configurations/'
%     folder as seperate files. The default case is 'YawCase3.m'.
%
%     - scriptOptions: this struct contains all simulation settings, not
%     related to the wind farm itself (solution methodology, outputs, etc.)
%
%     - OutputData (*): this struct contains all interesting outputs. It packs
%       several substructs, namely:
%       - *.Wp: this struct contains all the simulation settings related to 
%         the wind farm, the turbine inputs, the atmospheric properties, etc.
%         See WFSim.m for a more elaborate description of 'Wp'.
%
%       - *.sol_array: this cell array contains the system states at every
%         simulated time instant. Each cell entry contains a sol struct.
%         See WFSim.m for a more elaborate description of 'sol'. In
%         addition to the entries described in WFSim.m, each 'sol' struct
%         contains in addition:
%           *.sol.score: a struct containing estimation performance scores
%           such as the maximum estimation error, the RMS error, and the
%           computational cost (CPU time).
%           *.sol.measuredData: a struct containing the true (to be
%           estimated) values, and the measurement data fed into the
%           estimation algorithm.
%
%       - *.strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.
%
%%   Debugging and contributing:
%     - First, try to locate any errors by turning all possible outputs
%       on (printProgress, printConvergence, Animate, plotMesh).
%     - If you cannot solve your problems, reach out on the Github.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parallel Pool
myCluster = parcluster('local');
N_pool = myCluster.NumWorkers;
if isempty(gcp('nocreate')) ~= 0
    parpool(N_pool);
end

%% Define script settings
% Command window reporting settings
scriptOptions.printProgress     = 1;  % Print progress every timestep
scriptOptions.printConvergence  = 0;  % Print convergence parameters every timestep

% Visualization settings
scriptOptions.plotMesh          = 0;  % Show meshing and turbine locations
scriptOptions.Animate           = 10;  % Show results every x iterations (0: no plots)
   scriptOptions.plotContour    = 1;  % Show flow fields
   scriptOptions.plotPower      = 1;  % Plot true and predicted power capture vs. time
    scriptOptions.powerForecast = 0;  % Plot power forecast (0 = disabled, x = number of steps) (only if plotPower = 1)
   scriptOptions.plotError      = 0;  % plot RMS and maximum error vs. time
   scriptOptions.plotCenterline = 1;  % Plot centerline speed of the wake (m/s)

% Saving settings
scriptOptions.savePlots         = 1;  % Save all plots in external files at each time step
scriptOptions.saveWorkspace     = 1;  % Save complete workspace at the end of simulation
% scriptOptions.savePath          = ['/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs/results/tmp']; % Destination folder of saved files
scriptOptions.savePath          = ['/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs/results/tmp']; % Destination folder of saved files

% Configuration file
% configName = 'apc_9turb_alm_turb_dexkf_IFAC1DZ';
% configName = 'apc_9turb_adm_noturb';
configName = 'axi_2turb_alm_turb';
% axi_2turb_alm_turb_dexkf_IFAC_1DZ = axi_2turb_alm_turb + dexkf + fusion: IFAC
% + subsystem_length: 1D + typeCZ: Z
%% Execute the WFObs core code (+ overwrite meshing.m settings, if applicable)
WpOverwrite = struct(); % Struct to overwrite settings from meshing.m
WpOverwrite.sim.NN = 500; % Stop after [x] steps
run('WFObs_addpaths.m'); % Import libraries for WFObs & WFSim
% outputData = WFObs_core(scriptOptions,configName,WpOverwrite);
outputData = d_WFObs_core(scriptOptions,configName,WpOverwrite);