% Source files
scriptOptions.outputFilename = '2turb_adm_noturb';
scriptOptions.plotFrequency  = 1e9;         % Plot mapping every * instances (will always plot k == 1, set to high value for no plots after k == 1)
scriptOptions.sourcePath     = '/tudelft/ls/staff-group/3me/dcsc/DataDriven/Data/PALM/2turb_adm_matlab/noturb'

% Turbine properties directly from PALM or SOWFA. The reference frame is 
%   x (vertical, upwards pos.) - y (horizontal, rightwards pos.).
rawTurbData           = struct('Crx',[5700.0, 6456.0],'Cry',[1175.0, 1175.0]);
rawTurbData.Drotor    = [126.4, 126.4]; % Rotor diameter in (m)
rawTurbData.hubHeight = 90.0;           % Hub height in (m)

% Filtering
filterSettings.turbData.MM = true; % Apply moving average to turbine data
filterSettings.turbData.tL = 1;    % Window width to the left (seconds)
filterSettings.turbData.tR = 1;    % Window width to the right (seconds)

%filterSettings.Ur.nPts = 50;    % Number of rotor points (only if CT not available)
%filterSettings.Ur.MM   = true;  % Moving-average for rotor velocity (Only when CT not given)
%filterSettings.Ur.tL   = 5;     % Moving-average for rotor velocity (Only when CT not given)
%filterSettings.Ur.tR   = 5;     % Moving-average for rotor velocity (Only when CT not given)

filterSettings.CTp.MM   = true;  % Additional moving-mean average for Ct_prime
filterSettings.CTp.tL   = 3;     % Additional moving-mean average for Ct_prime
filterSettings.CTp.tR   = 3      % Additional moving-mean average for Ct_prime

% Desired output settings
meshSetup.dt          = 1.0 ; % Timestep in seconds
meshSetup.rho         = 1.20; % Air density (kg m^-3)
meshSetup.distance_S  = 180 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 750;  % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 250 ; % distance (m) west most left  turbine (from hub) to export
meshSetup.distance_E  = 250 ; % distance (m) east most right turbine (from hub) to export
meshSetup.Nx          = 40;   % Number of grid points in x-direction (-)
meshSetup.Ny          = 20    % Number of grid points in y-direction (-)