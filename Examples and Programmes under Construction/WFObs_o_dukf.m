function [Wp,sol,strucObs] = WFObs_o_dukf( strucObs,Wp,sys,sol,options) 

% Model Decomposition
tur = Wp.turbine.N;
if strucObs.stateEst || strucObs.measFlow
    stateLocArray = zeros(strucObs.size_output,2);
    for iii = 1:strucObs.size_output
        [~,loci,~]           = WFObs_s_sensors_nr2grid(iii,Wp.mesh);
        stateLocArray(iii,:) = [loci.x, loci.y];
    end
end
turbLocArray = zeros(Wp.turbine.N,2);
for iii = 1:Wp.turbine.N
    turbLocArray(iii,:) = [Wp.turbine.Crx(iii),Wp.turbine.Cry(iii)];
end
Subsys_length = strucObs.Subsys_length;
[x,d, tur,n, x_est,x_unest] = subsystem_turbine_nl(tur,stateLocArray,turbLocArray,Subsys_length);
    
%% Initialization step of the Unscented KF (at k == 1)

if sol.k==1 
    % Check UKF setting
    if strucObs.stateEst == 0 && strucObs.tune.est == 0
        error(['Please turn on state and/or parameter estimation. '...
            'Alternatively, select "sim" for open-loop simulations.']);
    end
    
    % Initialize state vector
    sol.x = [vec(sol.u(3:end-1,2:end-1)'); vec(sol.v(2:end-1,3:end-1)')];
    if options.exportPressures == 1 % Optional: add pressure terms
        sol.x = [sol.x; vec(sol.p(2:end-1,2:end-1)')];
        sol.x = sol.x(1:end-2); % Correction for how pressure is formatted
    end    
    
    % If do state estimation, then load covariance matrices
    if strucObs.stateEst
        x0 = sol.x;
        P0 = blkdiag(eye(Wp.Nu)*strucObs.P_0.u,eye(Wp.Nv)*strucObs.P_0.v);
        Qx = blkdiag(eye(Wp.Nu)*strucObs.Q_k.u,eye(Wp.Nv)*strucObs.Q_k.v);
        if options.exportPressures == 1 % Optional: add pressure terms
            P0 = blkdiag(Px,eye(Wp.Np)*strucObs.P_0.p);
            Qx = blkdiag(Qx,eye(Wp.Np)*strucObs.Q_k.p);
        end
    else
        % No state estimation, parameter only
        x0    = [];
        P0    = [];
        Qx    = [];
    end;
    
    % Add model parameters as states for online model adaption
%     if strucObs.tune.est
%         for iT = 1:length(strucObs.tune.vars)
%             tuneP       = strucObs.tune.vars{iT};
%             dotLoc      = findstr(tuneP,'.');
%             subStruct   = tuneP(1:dotLoc-1);
%             structVar   = tuneP(dotLoc+1:end);
%             x0          = [x0; Wp.(subStruct).(structVar)];
%             P0          = blkdiag(P0,strucObs.tune.P_0(iT));
%             Qx          = blkdiag(Qx,strucObs.tune.Q_k(iT));
% 
%             % Save to strucObs for later
%             strucObs.tune.subStruct{iT} = subStruct;
%             strucObs.tune.structVar{iT} = structVar;
%         end
%     end
    
    % Calculate initial ensemble
    L               = length(x0);
    lambda          = strucObs.alpha^2*(L+strucObs.kappa)-L;
    gamma           = sqrt(L+lambda);
    
    strucObs.nrens  = 2*L+1; % number of sigma points
    strucObs.nrobs  = length(strucObs.obs_array); % number of measurements  
    
    % Weight vectors for UKF update
    strucObs.Wm     = [lambda/(L+lambda); repmat(1/(2*(L+lambda)),2*L,1)];
    strucObs.Wc     = [lambda/(L+lambda)+(1-strucObs.alpha^2+strucObs.beta); ...
                       repmat(1/(2*(L+lambda)),2*L,1)];
    tempMat         = eye(2*L+1)-repmat(strucObs.Wm,1,2*L+1);
    strucObs.W      = tempMat*diag(strucObs.Wc)*tempMat';
               
    % Covariance matrices for UKF
    strucObs.Qx = Qx;
    strucObs.Pk = P0;
    strucObs.Rx = eye(strucObs.nrobs)*strucObs.R_k;
    if strucObs.measPw 
        strucObs.Rx = blkdiag(strucObs.Rx,eye(Wp.turbine.N)*strucObs.R_ePw);
    end
    
    % Other UKF-related parameters
    strucObs.L      = L;
    strucObs.lambda = lambda;
    strucObs.gamma  = gamma;
    strucObs.M      = strucObs.nrobs+strucObs.measPw*Wp.turbine.N; % total length of measurements
end

%%%%%%%%%%%%
% Recursive Steps

% Calculate sigma points
if strucObs.stateEst % Write the system states to the sigma points
    strucObs.Aen = repmat(sol.x,1,strucObs.nrens); 
else
    strucObs.Aen = [];
end;

% Append the sigma points with model parameters
% if strucObs.tune.est
%     for iT = 1:length(strucObs.tune.vars) 
%         strucObs.Aen = [strucObs.Aen; repmat(Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}),1,strucObs.nrens)];
%     end
% end

% Distribute sigma points around the mean
Uscented_devs                    = strucObs.gamma*sqrt(strucObs.Pk);
strucObs.Aen(:,2:strucObs.L+1)   = strucObs.Aen(:,2:strucObs.L+1)   + Uscented_devs;
strucObs.Aen(:,strucObs.L+2:end) = strucObs.Aen(:,strucObs.L+2:end) - Uscented_devs;

%% Parallelized solving of the forward propagation step in the UKF
% tic
Aenf    = cell(tur,1);
Yenf    = cell(tur,1);
L       = strucObs.L;
obs_array = cell(tur,1);
for i = 1:tur
    nrens{i} = [1,1+x{i}',(1+L)+x{i}'];
    Aenf{i}  = zeros(length(x{i}),length(nrens{i}));   % Initialize empty forecast matrix
    Yenf{i}  = zeros(strucObs.M,length(nrens{i}));   % Initialize empty output matrix
    Atmp = Aenf{i};
    Ytmp = Yenf{i};
    for j = 1:length(x{i})
        if sum(double(x{i}(j) == strucObs.obs_array)) == 1
            obs_array{i} = [obs_array{i},x{i}(j)];
        end
    end
    parfor(ji=1:nrens{i})
        syspar   = sys; % Copy system matrices
        solpar   = sol; % Copy solution from prev. time instant
        Wppar    = Wp;  % Copy meshing struct

        % Import solution from sigma point
        if strucObs.stateEst
            % Load sigma point as solpar.x
            solpar.x   = strucObs.Aen(1:strucObs.size_output,ji);
    %         solpar.x   = strucObs.Aen(x{i},ji);
            [solpar,~] = MapSolution(Wppar,solpar,Inf,options);
        end;

        % Update Wp with values from the sigma points
    %     if strucObs.tune.est
    %         tuneParam_tmp = zeros(length(strucObs.tune.vars),1);
    %         for iT = 1:length(strucObs.tune.vars)
    %             % Threshold using min-max to avoid crossing lb/ub
    %             tuneParam_tmp(iT) = min(strucObs.tune.ub(iT),max(strucObs.tune.lb(iT),...
    %                                 strucObs.Aen(end-length(strucObs.tune.vars)+iT,ji)));
    %             Wppar.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = tuneParam_tmp(iT);
    %         end
    %     end

        % Forward propagation
        solpar.k          = solpar.k - 1;
        [ solpar,syspar ] = WFSim_timestepping( solpar, sys, Wppar, options );

        xf = [];
        if strucObs.stateEst 
            xf = [xf; solpar.x(x{i},1)];
        end
    %     if strucObs.tune.est
    %         xf = [xf; tuneParam_tmp];
    %     end

        % Write forecasted state to ensemble forecast matrix   
        Atmp(:,ji) = xf;

        % Calculate output vector
        if strucObs.measPw
            Ytmp(:,ji) = [solpar.x(strucObs.obs_array); solpar.turbine.power];
        else
            Ytmp(:,ji) = [solpar.x(strucObs.obs_array)];
        end    
    end
    Aenf{i} = Atmp;
    Yenf{i} = Ytmp;
end

% parfor(ji=1:strucObs.nrens)
%     syspar   = sys; % Copy system matrices
%     solpar   = sol; % Copy solution from prev. time instant
%     Wppar    = Wp;  % Copy meshing struct
%    
%     % Import solution from sigma point
%     if strucObs.stateEst
%         % Load sigma point as solpar.x
%         solpar.x   = strucObs.Aen(1:strucObs.size_output,ji);
%         [solpar,~] = MapSolution(Wppar,solpar,Inf,options);
%     end
% 
%     % Forward propagation
%     solpar.k          = solpar.k - 1;
%     [ solpar,syspar ] = WFSim_timestepping( solpar, sys, Wppar, options );
%     
%     xf = [];
%     if strucObs.stateEst 
%         xf = [xf; solpar.x(1:strucObs.size_output,1)];
%     end
% 
%     % Write forecasted state to ensemble forecast matrix   
%     Aenf(:,ji) = xf;
%     
%     % Calculate output vector
%     if strucObs.measPw
%         Yenf(:,ji) = [solpar.x(strucObs.obs_array); solpar.turbine.power];
%     else
%         Yenf(:,ji) = [solpar.x(strucObs.obs_array)];
%     end    
% end
% toc

%% Analysis update of the Unscented KF
% tic
if strucObs.measPw
    y_meas = [sol.measuredData.sol(strucObs.obs_array);sol.measuredData.power];
else
    y_meas = [sol.measuredData.sol(strucObs.obs_array)];
end
% xmean   = cell(tur,1);
% ymean   = cell(tur,1);
% Aenft   = cell(tur,1);
% Yenft   = cell(tur,1);
z       = cell(tur,1);
Z       = cell(tur,1);
parfor i=1:tur
    xmean = sum(repmat(strucObs.Wm(nrens{i},1)',length(x{i}),1) .* Aenf{i}, 2);
    ymean = sum(repmat(strucObs.Wm(nrens{i},1)',strucObs.M,1) .* Yenf{i}, 2);
    Aenft = Aenf{i}-repmat(xmean,1,length(nrens{i}));
    Yenft = Yenf{i}-repmat(ymean,1,length(nrens{i}));
    Pfxxk = Aenft*strucObs.W(nrens{i},nrens{i})*Aenft' + strucObs.Qx(x{i},x{i}); % Pxx for k|k-1
    Pfyyk = Yenft*strucObs.W(nrens{i},nrens{i})*Yenft' + strucObs.Rx; % Pxy for k|k-1
    Pfxyk = Aenft*strucObs.W(nrens{i},nrens{i})*Yenft';               % Pyy for k|k-1

    Kk      = Pfxyk * pinv(Pfyyk);
    z{i}    = xmean + Kk*(y_meas-ymean);
    Z{i}    = Pfxxk - Kk * Pfyyk * Kk';
end
% toc
typeCZ = 1;       % 1 if Z = Co-Variance, 2 if Z = Information, 
[xtmp,Ptmp] = fuze(z,Z,x,tur,L,x_est,typeCZ);
z{i};
Z{i};

xSolAll             = zeros(L,1);
xSolAll(x_est,1)    = xtmp;
xSolAll(x_unest,1)  = sol.x(x_unest);

strucObs.Px                 = 2*eye(L,L);
strucObs.Px(x_est,x_est)    = Ptmp;
%     xSolAll     = xmean + Kk*(y_meas-ymean);
%     strucObs.Px = Pfxxk - Kk * Pfyyk * Kk';

%% Post-processing
% if strucObs.tune.est
%     % Update model parameters with the optimal estimate
%     for iT = 1:length(strucObs.tune.vars) % Write optimally estimated values to Wp
%         Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}) = min(...
%             strucObs.tune.ub(iT),max(strucObs.tune.lb(iT),xSolAll(end-length(strucObs.tune.vars)+iT)));
%     end
% end

% Update states, either from estimation or through open-loop
if strucObs.stateEst
    sol.x    = xSolAll(1:strucObs.size_output); % Write optimal estimate to sol
    [sol,~]  = MapSolution(Wp,sol,Inf,options); % Map solution to flow fields
    [~,sol]  = Actuator(Wp,sol,options);        % Recalculate power after analysis update
else
    % Note: this is identical to 'sim' case in WFObs_o(..)
    sol.k    = sol.k - 1; % Necessary since WFSim_timestepping(...) already includes time propagation
    [sol,~]  = WFSim_timestepping(sol,sys,Wp,options);
end

