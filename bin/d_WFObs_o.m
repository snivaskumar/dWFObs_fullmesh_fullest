function [ Wp,sol,strucObs ] = d_WFObs_o(strucObs,Wp,sys,sol,options)
% WFOBS_O  Header function to call the correct estimation function
%
%   SUMMARY
%    This code calls the correct estimation function (EnKF, UKF, ExKF, ...)
%    based on what is specified in strucObs.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.
%
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%         Wp.Nu:      Number of model states concerning longitudinal flow.
%         Wp.Nv:      Number of model states concerning lateral flow.
%         Wp.Np:      Number of model states concerning pressure terms.
%         Wp.sim:     Substruct containing timestep and simulation length.
%         Wp.turbine: Substruct containing turbine properties and settings.
%         Wp.site:    Substruct containing freestream atmospheric properties.
%         Wp.mesh:    Substruct containing topology and meshing settings.
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.A:     System matrix A in the grand picture: A*sol.x = b
%         sys.b:     System vector b in the grand picture: A*sol.x = b
%         sys.pRCM:  Reverse Cuthill-McKee algorithm for solving A*x=b faster.
%         sys.B1:    Important matrix in the boundary conditions.
%         sys.B2:    Important matrix in the boundary conditions.
%         sys.bc:    Important vector in the boundary conditions.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.k:     Discrete timestep  to which these system states belong
%         sol.time:  Actual time (in s) to which these system states belong
%         sol.x:     True system state (basically flow field excluding bcs)
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.p:     Instantaneous pressure field over the mesh (in Pa)
%         sol.uu:    Same as sol.u, used for convergence
%         sol.vv:    Same as sol.v, used for convergence
%         sol.pp:    Same as sol.p, used for convergence
%         sol.turbine: a struct containing relevant turbine outputs such as
%         the ax. ind. factor, the generated power, and the ct coefficient
%         sol.measuredData: a struct containing the true flow field and the
%         measurements used for estimation
%         sol.score: a struct containing estimator performance measures
%         such as the magnitude and location of the maximum flow estimation
%         error, the RMSE, and the computational cost (CPU time)
%
%     - options: this struct contains all simulation settings, not related
%       to the wind farm itself (solution methodology, outputs, etc.). See
%       also 'scriptOptions' (options == scriptOptions).
%       
switch lower(strucObs.filtertype)
    case 'dexkf'
        % Distributed Extended Kalman filtering
        [Wp,sol,strucObs] = d_WFObs_o_dexkf(strucObs,Wp,sys,sol,options);
    case 'exkf'
        % Extended Kalman filtering
        [Wp,sol,strucObs] = d_WFObs_o_exkf(strucObs,Wp,sys,sol,options);
    case 'smo'
        % Sliding mode observer
        [Wp,sol,strucObs] = WFObs_o_smo(strucObs,Wp,sys,sol,options);        
    case 'enkf'
        % Ensemble Kalman filtering
        [Wp,sol,strucObs] = WFObs_o_enkf(strucObs,Wp,sys,sol,options);
    case 'ukf'
        % Unscented Kalman filtering
%         tic
        [Wp,sol,strucObs] = WFObs_o_ukf( strucObs,Wp,sys,sol,options); 
%         toc
    case 'dukf'
        % Unscented Kalman filtering
%         tic
        [Wp,sol,strucObs] = WFObs_o_dukf( strucObs,Wp,sys,sol,options); 
%         toc
    case 'sim'
        % Open-loop simulations
        sol.k    = sol.k - 1; % Necessary since WFSim_timestepping(...) already includes time propagation
        sol.time = Wp.sim.time(sol.k+1);% Timestep backward
        [sol,~]  = d_WFSim_timestepping(sol,sys,Wp,options);
    otherwise
        error('not a valid filter specified.');
end;

end