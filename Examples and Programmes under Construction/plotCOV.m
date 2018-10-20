function plotCOV( Wp,sol_array,sys,scriptOptions,strucObs);

sol             = sol_array(end);
measuredData    = sol.measuredData;
% C               = strucObs.Pk;
Nx              = Wp.mesh.Nx;
Ny              = Wp.mesh.Ny;
% 
x = [vec(sol.u(3:end-1,2:end-1)'); ...
     vec(sol.v(2:end-1,3:end-1)')];
y = [vec(measuredData.uq(3:end-1,2:end-1)'); ...
     vec(measuredData.vq(2:end-1,3:end-1)')];
e = x - y; 
% % MDis    = e'*inv(C)*e; 
% MDis    = mahal(y,x); 

% ex_Wp = Wp;
% ex_sol_array = sol_array;
% ex_sys = sys;
% ex_scriptOptions = scriptOptions;
% ex_strucObs = strucObs;
% 
% address;
% 
% Wp = ex_Wp;
% sol_array = ex_sol_array;
% sys = ex_sys;
% scriptOptions = ex_scriptOptions;
% strucObs = ex_strucObs;


Nx  = Wp.mesh.Nx;
Ny  = Wp.mesh.Ny;
if ~length(strfind(scriptOptions.savePath,'EnKF'))
    PP          = strucObs.Pk;
    P           = diag(strucObs.Pk);
else
    Aenf    = strucObs.Aen;
    Aenft   = Aenf-repmat(mean(Aenf,2),1,strucObs.nrens);
%     P       = ( 1/(strucObs.nrens-1) )*Aenft*Aenft';
    PP      = Aenft*Aenft';
    P       = diag(PP);
end
tmp                     = max(max(P));
Puv                     = tmp*ones(Nx,Ny);
Puv(3:end-1,2:end-1)    = reshape(P(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
[Puv,~,~] = Updateboundaries(Nx,Ny,Puv,Puv,Puv);    

% % Ptmp = min(min(Puv));
% Ptmp = mean(mean(Puv));
% for i = 1:length(Ix)
%     Puv(Ix(i),Iy(i)) = Ptmp;
% end

C       = P;
error   = e.*e;
MDis    = error./C; 
tmp                     = max(max(MDis));
dis                     = tmp*ones(Nx,Ny);
dis(3:end-1,2:end-1)    = reshape(MDis(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
[dis,~,~] = Updateboundaries(Nx,Ny,dis,dis,dis);

hFigs = {};
scrsz = get(0,'ScreenSize'); 
hFigs{1}=figure('color',[1 1 1],'Position',[50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)], 'MenuBar','none','ToolBar','none','visible', 'on');
set(hFigs{1},'defaultTextInterpreter','latex')
data{1} = struct('x',Wp.mesh.ldyy, 'y',Wp.mesh.ldxx2,'z',Puv,'title','Co_variance');
data{2} = struct('x',Wp.mesh.ldyy, 'y',Wp.mesh.ldxx2,'z',dis,'title','Mahalanobis Distance');
        
% applied correction for yaw angle: wake was forming at wrong side
rotorRotation = -.5*Wp.turbine.Drotor*exp(1i*-Wp.turbine.input(sol.k).phi'*pi/180); 
    
% Plot velocities in a contourf figure
set(0,'CurrentFigure',hFigs{1}); clf
subplotDim = numSubplots(length(data));
for j = 1:length(data)
	subplot(subplotDim(1),subplotDim(2),j);
    if j == 1
%         cmax = 1;
        cmax = (max(max(Puv)));
        cmin = (min(min(Puv)));
    end
    if j == 2
%         cmax = (max(max(dis)));
        cmax = 3;
        cmin = (min(min(dis)));
    end
    contourf(data{j}.x,data{j}.y,data{j}.z,cmin:0.1:cmax,'Linecolor','none');
    title([data{j}.title ' (t = ' num2str(sol.time) ')'])
    hold all; colorbar;
    caxis([cmin cmax]);
    axis equal; axis tight;
    xlabel('y-direction')
    ylabel('x-direction')         
    hold all
    % Turbines
    for kk=1:Wp.turbine.N
        Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
        Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
                rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                'FaceColor','w')                  
        plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
        plot(Qy,Qx,'w','linewidth',2)              
    end
    % Sensors
    if strucObs.measFlow
    	plot([strucObs.obs_array_locu.y],[strucObs.obs_array_locu.x],'wo','lineWidth',3.0,'displayName','Sensors');
        plot([strucObs.obs_array_locu.y],[strucObs.obs_array_locu.x],'ro','displayName','Sensors');
    end
	set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
end;
colormap(jet)
end