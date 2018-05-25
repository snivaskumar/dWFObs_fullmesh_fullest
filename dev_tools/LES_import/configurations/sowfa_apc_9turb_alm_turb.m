% Source files
scriptOptions.outputFilename = 'apc_9turb_alm_turb';
scriptOptions.plotFrequency  = 1e9;   % Plot mapping every * instances (will always plot k == 1, set to high value for no plots after k == 1)
scriptOptions.sourcePath     = '/tudelft/ls/staff-group/3me/dcsc/DataDriven/Data/SOWFA/APC_12_wake_noyaw_pitch_agc_50/sliceDataInstant'

% Turbine properties directly from PALM or SOWFA. The reference frame is 
%   x (vertical, upwards pos.) - y (horizontal, rightwards pos.).
tLocs = [1465.9, 898.1  ; ... % turbine 1
         1276.3, 1226.5 ; ... % turbine 2
         1086.7, 1554.9 ; ... % turbine 3
         2013.2, 1214.1 ; ... % turbine 4
         1823.6, 1542.5 ; ... % turbine 5
         1634.0, 1870.9 ; ... % turbine 6
         2560.5, 1530.1 ; ... % turbine 7
         2370.9, 1858.5 ; ... % turbine 8
         2181.3, 2186.9 ];    % turbine 9
rawTurbData           = struct('Crx',tLocs(:,2),'Cry',tLocs(:,1));			   
rawTurbData.Drotor    = 126.4*ones(1,9); % Rotor diameter in (m)
rawTurbData.hubHeight = 90.0             % Hub height in (m)

% Filtering
filterSettings.turbData.MM = true; % Apply moving average to turbine data
filterSettings.turbData.tL = 1;    % Window width to the left (seconds)
filterSettings.turbData.tR = 1;    % Window width to the right (seconds)

filterSettings.Ur.nPts = 50;    % Number of rotor points (only if CT not available)
filterSettings.Ur.MM   = true;  % Moving-average for rotor velocity (Only when CT not given)
filterSettings.Ur.tL   = 5;     % Moving-average for rotor velocity (Only when CT not given)
filterSettings.Ur.tR   = 5;     % Moving-average for rotor velocity (Only when CT not given)

filterSettings.CTp.MM   = true;  % Additional moving-mean average for Ct_prime
filterSettings.CTp.tL   = 3;     % Additional moving-mean average for Ct_prime
filterSettings.CTp.tR   = 3      % Additional moving-mean average for Ct_prime

% Desired output settings
meshSetup.dt          = 1.0 ; % Timestep in seconds
meshSetup.rho         = 1.20; % Air density (kg m^-3)
meshSetup.distance_S  = 400 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 850 ; % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 400 ; % distance (m) west most left  turbine (from hub) to export
meshSetup.distance_E  = 400 ; % distance (m) east most right turbine (from hub) to export
meshSetup.Nx          = 100 ; % Number of grid points in x-direction (-)
meshSetup.Ny          = 42   % Number of grid points in y-direction (-)