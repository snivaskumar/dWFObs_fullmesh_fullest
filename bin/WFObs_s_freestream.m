function [ Wp,sol,sys,strucObs ] = WFObs_s_freestream( Wp,sol,sys,strucObs )
% WFOBS_S_FREESTREAM  Estimate the freestream conditions u_Inf and v_Inf
%
%   SUMMARY
%    This code estimates the freestream flow speed and direction in 2D
%    using the power measurements from the most upstream row of turbines.
%    It uses wind vane measurements to determine the WD and thereby the
%    most upstream turbines. Then, it uses power to determine WS.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%         Wp.site:    Substruct containing freestream atmospheric properties.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.uu:    Same as sol.u, used for convergence
%         sol.vv:    Same as sol.v, used for convergence
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.B1:    Important matrix in the boundary conditions.
%         sys.B2:    Important matrix in the boundary conditions.
%         sys.bc:    Important vector in the boundary conditions.
%
%       - *.strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.
%
if strucObs.U_Inf.estimate && sol.k > 25 % Skip initialization period   
    % Import variables
    input        = Wp.turbine.input(sol.k);
    measuredData = sol.measuredData;
    
    wd = 270.; % wind direction in degrees. should actually be something like:
    % wd = mean(measured.windVaneMeasurements)
    % but currently no anemometer measurements yet from SOWFA...
    
    upstreamTurbines = WFObs_s_freestream_findUnwaked( Wp, wd );
    
    U_est = [];
    eta = 0.95; % Correction factor to correct for ADM/WFSim mismatch
    psc = Wp.turbine.powerscale; % Powerscale [-]
    Rho = Wp.site.Rho; % Air density [kg/m3]
    Ar  = 0.25*pi*Wp.turbine.Drotor^2; % Rotor swept area [m2]
    CTp  = [Wp.turbine.input(sol.k).CT_prime]; % CT' [-]

    % Calculate U_Inf for each turbine
    U_Inf_vec = (1+0.25*CTp(upstreamTurbines)).*((measuredData.power(upstreamTurbines)...
                ./(eta*psc*0.5*Rho*Ar*CTp(upstreamTurbines))).^(1/3));
    
    % Determine previous average and current instantaneous U_Inf
    U_Inf_Previous      = sqrt(sol.u(1,1)^2+sol.v(1,1)^2); % = sqrt(Wp.site.u_Inf^2+Wp.site.v_Inf^2);
    U_Inf_Instantaneous = mean(U_Inf_vec);
    
    % Low-pass filter the mean instantaneous U_Inf
    tau = strucObs.U_Inf.intFactor; % Time constant of LPF. 0 = instant updates, 1 = no updates. 0.86 means after 30 seconds, 1% of old solution is left
    U_Inf_Filtered = tau*U_Inf_Previous+(1-tau)*U_Inf_Instantaneous;
    
    % Reformat to x and y direction
    u_Inf = U_Inf_Filtered*cosd(270-wd);
    v_Inf = U_Inf_Filtered*sind(270-wd);
    
    % Shift entire solution space to accomodate for new freestream conditions
    [sol.u,sol.uu] = deal(sol.u+(u_Inf-Wp.site.u_Inf)); % Update all states
    [sol.v,sol.vv] = deal(sol.v+(v_Inf-Wp.site.v_Inf)); % Update all states

    % Update inflow parameters in Wp
    Wp.site.u_Inf = u_Inf;
    Wp.site.v_Inf = v_Inf;
        
    % Compute system boundary conditions and corresponding matrices B1, B2
    [sys.B1,sys.B2,sys.bc] = Compute_B1_B2_bc(Wp); 
end
end