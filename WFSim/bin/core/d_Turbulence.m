% Define spatial varying mixing-length parameter
% Should work for wd = 270deg
% if N == 1
%     if sum(ismember(ID,upstreamTurbines)) == 0
%         mixing_length    = lmu*0.5*Drotor*ones(Nx,Ny);
%     else
%         x                = [zeros(1,xline(1)+n) linspace(0,lmu,Nx-xline(1)-n)]';
%         y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
%         mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
%     end
% elseif N == 2
%     if sum(ismember(ID,upstreamTurbines)) == 0
%         mixing_length    = lmu*0.5*Drotor*ones(Nx,Ny);
%     elseif sum(ismember(ID,upstreamTurbines)) == 1     % 2 turbine case
%         x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(2)-xline(1)-4*n)]';
%         x                = [x;zeros(4*n,1);linspace(0,lmu,Nx-xline(2)-n)'];
%         y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
%         mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
%     elseif (sum(ismember(ID,upstreamTurbines)) == 2) 
%         xline            = sort(unique(xline));
%         x                = [zeros(1,xline(1)+n) linspace(0,lmu,Nx-xline(1)-n)]';
%         y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,yline{2}(1)-yline{1}(end)-1)...
%                          ones(1,length(yline{2})) zeros(1,Ny-yline{2}(end))] ;
%         mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
%     end
% elseif N == 3
%     if sum(ismember(ID,upstreamTurbines)) == 0
%         mixing_length    = lmu*0.5*Drotor*ones(Nx,Ny);
%     elseif sum(ismember(ID,upstreamTurbines)) == 1     % 3 turbine case
%         disp('Sub-system will not be this big.')
%     else
%         xline            = sort(unique(xline));
%         x                = [zeros(1,xline(1)+n) linspace(0,lmu,Nx-xline(1)-n)]';
%         y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,yline{2}(1)-yline{1}(end)-1) ...
%                         ones(1,length(yline{2})) zeros(1,yline{3}(1)-yline{2}(end)-1) ...
%                         ones(1,length(yline{3})) zeros(1,Ny-yline{3}(end))];
%         mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor; 
%     end
% elseif N == 4
%     if sum(ismember(ID,upstreamTurbines)) == 0
%         mixing_length    = lmu*0.5*Drotor*ones(Nx,Ny);
%     else % (if sum(ismember(ID,upstreamTurbines)) == 2)
%         xline            = sort(unique(xline));
%         x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(1+1)-xline(1)-4*n)]';
%         x                = [x;zeros(4*n,1);linspace(0,lmu,Nx-xline(2)-n)'];
%         y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,yline{2}(1)-yline{1}(end)-1)...
%                          ones(1,length(yline{2})) zeros(1,Ny-yline{2}(end))] ;
%         mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
%     end
% elseif N == 6
%     if sum(ismember(ID,upstreamTurbines)) == 0
%         mixing_length    = lmu*0.5*Drotor*ones(Nx,Ny);
%     else
%         xline            = sort(unique(xline));
%         x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(2)-xline(1)-m)]';
%         x                = [x;zeros(m,1);linspace(0,lmu,xline(3)-xline(2)-m)'];
%         x                = [x;zeros(m,1);linspace(0,lmu,Nx-xline(3)-n)'];
%         y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,yline{2}(1)-yline{1}(end)-1)...
%                              ones(1,length(yline{2})) zeros(1,Ny-yline{2}(end))] ;
%         mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
%     end
% elseif N == 9
%     xline            = sort(unique(xline));
%     x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(2)-xline(1)-m)]';
%     x                = [x;zeros(m,1);linspace(0,lmu,xline(3)-xline(2)-m)'];
%     x                = [x;zeros(m,1);linspace(0,lmu,Nx-xline(3)-n)'];
%     y                = [zeros(1,yline{3}(1)-1) ones(1,length(yline{3})) zeros(1,yline{5}(1)-yline{3}(end)-1) ...
%                         ones(1,length(yline{5})) zeros(1,yline{7}(1)-yline{5}(end)-1) ...
%                         ones(1,length(yline{7})) zeros(1,Ny-yline{7}(end))];
%     mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
% else
%     mixing_length    = lmu*0.5*Drotor*ones(Nx,Ny);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mixing Length matrix for the subsystems are created 
% by creating Mixing Length matrix for the entire mesh 
% & cropping the part for the subsystem by using the spatial co-ordinates
Turbulence = upper(Turbulence);
if strcmp(Turbulence,'CENTRALIZE')
    tur     = Wp.tur;
    NNx     = Wp.actualmesh.Nx;
    NNy     = Wp.actualmesh.Ny;
    Nxb     = Wp.mesh.Nxb;
    Nxe     = Wp.mesh.Nxe;
    Nyb     = Wp.mesh.Nyb;
    Nye     = Wp.mesh.Nye;
    NNxb    = cell2mat(Wp.mesh.NNNxb);
    NNxe    = cell2mat(Wp.mesh.NNNxe);
    NNyb    = cell2mat(Wp.mesh.NNNyb);
    NNye    = cell2mat(Wp.mesh.NNNye);
    xxline  = Wp.mesh.xxline;
    yyline  = Wp.mesh.yyline;
    xxline  = cell2mat(xxline);
    for i = 1:tur
        yyline{i} = cell2mat(yyline{i});
    end
    if tur==1
        x                = [zeros(1,xxline(1)+n) linspace(0,lmu,NNx-xxline(1)-n)]';
        y                = [zeros(1,yyline{1}(1)-1) ones(1,length(yyline{1})) zeros(1,NNy-yyline{1}(end))] ;
        mixing_length    = (repmat(x,1,NNy).*repmat(y,NNx,1))*0.5*Drotor;
    elseif tur==2
        x                = [zeros(1,xxline(1)+n) linspace(0,lmu,xxline(1+1)-xxline(1)-4*n)]';
        x                = [x;zeros(4*n,1);linspace(0,lmu,NNx-xxline(2)-n)'];
        y                = [zeros(1,yyline{1}(1)-1) ones(1,length(yyline{1})) zeros(1,NNy-yyline{1}(end))] ;
        mixing_length    = (repmat(x,1,NNy).*repmat(y,NNx,1))*0.5*Drotor;
    elseif tur==3
        xxline            = sort(unique(xxline));
        x                = [zeros(1,xxline(1)+n) linspace(0,lmu,xxline(2)-xxline(1)-m)]';
        x                = [x;zeros(m,1);linspace(0,lmu,xxline(3)-xxline(2)-m)'];
        x                = [x;zeros(m,1);linspace(0,lmu,NNx-xxline(3)-n)'];
        y                = [zeros(1,yyline{1}(1)-1) ones(1,length(yyline{1})) zeros(1,NNy-yyline{1}(end))] ;
        mixing_length    = (repmat(x,1,NNy).*repmat(y,NNx,1))*0.5*Drotor;
    elseif tur==6
        xxline            = sort(unique(xxline));
        x                = [zeros(1,xxline(1)+n) linspace(0,lmu,xxline(2)-xxline(1)-m)]';
        x                = [x;zeros(m,1);linspace(0,lmu,xxline(3)-xxline(2)-m)'];
        x                = [x;zeros(m,1);linspace(0,lmu,NNx-xxline(3)-n)'];
        y                = [zeros(1,yyline{1}(1)-1) ones(1,length(yyline{1})) zeros(1,yyline{2}(1)-yyline{1}(end)-1)...
                             ones(1,length(yyline{2})) zeros(1,NNy-yyline{2}(end))] ;
        mixing_length    = (repmat(x,1,NNy).*repmat(y,NNx,1))*0.5*Drotor;
    elseif tur==9
        xxline            = sort(unique(xxline));
        x                = [zeros(1,xxline(1)+n) linspace(0,lmu,xxline(2)-xxline(1)-m)]';
        x                = [x;zeros(m,1);linspace(0,lmu,xxline(3)-xxline(2)-m)'];
        x                = [x;zeros(m,1);linspace(0,lmu,NNx-xxline(3)-n)'];
        y                = [zeros(1,yyline{3}(1)-1) ones(1,length(yyline{3})) zeros(1,yyline{5}(1)-yyline{3}(end)-1) ...
                            ones(1,length(yyline{5})) zeros(1,yyline{7}(1)-yyline{5}(end)-1) ...
                            ones(1,length(yyline{7})) zeros(1,NNy-yyline{7}(end))];
        mixing_length    = (repmat(x,1,NNy).*repmat(y,NNx,1))*0.5*Drotor;
    else
        mixing_length    = lmu*0.5*Drotor*ones(NNx,NNy);
    end 
    mixing_length = mixing_length(Nxb:Nxe,Nyb:Nye);  
    ab = zeros(size(mixing_length));
    Ix = find(xxline==xline);
    for i = 1:tur
        if length(yyline{i})==length(yline{1})
            if sum( yyline{i}==yline{1} ) == length( yline{1} )
                Iy = i;
            end
        end
    end
    if Ix < length(xxline)
        ab(xline:xxline(Ix+1)-1,yyline{Iy}) = mixing_length(xline:xxline(Ix+1)-1,yyline{Iy});
    else
        ab(xline:end,yyline{Iy}) = mixing_length(xline:end,yyline{Iy});
    end
    mixing_length = ab;
else
    if N==1
        x                = [zeros(1,xline(1)+n) linspace(0,lmu,Nx-xline(1)-n)]';
        y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
        mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
    else
        mixing_length    = lmu*0.5*Drotor*ones(Nx,Ny);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(mixing_length,1)>1
    H                = fspecial('disk',1); % You need Nx,Nx to keep it symmetric
    mixing_length    = filter2(H,mixing_length);
end

switch lower(Wp.site.Turbulencemodel)
    
    case lower('WFSim5')
        
        % For u-momentum equation
        ax.Tex              = zeros(Nx,Ny);
        ax.Twx              = zeros(Nx,Ny);
        ax.Tnx              = zeros(Nx,Ny);
        ax.Tsx              = zeros(Nx,Ny);
        
        ax.Tnex             = zeros(Nx,Ny);
        ax.Tsex             = zeros(Nx,Ny);
        
        %1
        ax.Tex(1:Nx-1,2:Ny-1) = Rho*(mixing_length(1:Nx-1,2:Ny-1).^2).*(dyy2(1:Nx-1,2:Ny-1)./(dyy(1:Nx-1,2:Ny-1).*dxx2(1:Nx-1,2:Ny-1)))...
            .*abs(u(1:Nx-1,2:Ny-1)-u(1:Nx-1,1:Ny-2));
        ax.Twx(2:Nx,2:Ny-1)   = Rho*(mixing_length(2:Nx,2:Ny-1).^2).*(dyy2(1:Nx-1,2:Ny-1)./(dyy(2:Nx,2:Ny-1).*dxx2(2:Nx,2:Ny-1)))...
            .*abs(u(1:Nx-1,1:Ny-2)-u(1:Nx-1,2:Ny-1));
        
        ax.aE             = ax.aE + ax.Tex;
        ax.aW             = ax.aW + ax.Twx;
        ax.aP             = ax.aP + ax.Twx + ax.Tex;
        
        %2
        ax.Tnx(1:Nx,1:Ny-1)= Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,1:Ny-1).^2))...
            .*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1));
        ax.Tsx(1:Nx,2:Ny)  = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2))...
            .*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny));
        
        ax.aN             = ax.aN + ax.Tnx;
        ax.aS             = ax.aS + ax.Tsx;
        ax.aP             = ax.aP + ax.Tnx + ax.Tsx;
        
        ax.Tnex1(1:Nx,1:Ny-1)= .5*Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(1./(dyy(1:Nx,1:Ny-1)))...
            .*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1))...
            +Wp.site.u_Inf*dxx(1:Nx,2:Ny)/2;
        ax.Tsex1(1:Nx,2:Ny)  = .5*Rho*(mixing_length(1:Nx,2:Ny).^2).*(1./(dyy(1:Nx,2:Ny)))...
            .*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny))...
            +Wp.site.u_Inf*dxx(1:Nx,2:Ny)/2;
        
        ax.Tnex2(1:Nx,1:Ny-1)= .5*Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(1./(dyy(1:Nx,1:Ny-1)))...
            .*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1))...
            +Wp.site.u_Inf*dxx(1:Nx,2:Ny)/2;
        ax.Tsex2(1:Nx,2:Ny)  = .5*Rho*(mixing_length(1:Nx,2:Ny).^2).*(1./(dyy(1:Nx,2:Ny)))...
            .*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny))...
            +Wp.site.u_Inf*dxx(1:Nx,2:Ny)/2;
        
        Axo  = sparse((Ny-2)*(Nx-3),(Nx-2)*(Ny-3));
        Bx_W = sparse(Ny-2,Ny-3); Bx_E   = sparse(Ny-2,Ny-3);
        
        for x= 3:Nx-1;
            swx      = spdiags(ax.Tsex1(x,3:Ny-1)',0,Ny-3,Ny-3);%-Wp.site.u_Inf*dxx(1:Ny-3,1:Ny-3)/2;%v_{I-1,j}
            nwx      = -spdiags(ax.Tnex1(x,2:Ny-2)',0,Ny-3,Ny-3);%+Wp.site.u_Inf*dxx(1:Ny-3,1:Ny-3)/2;    %v_{I-1,j+1}
            
            Bx_W     = [sparse(1,Ny-3); swx] + [nwx; sparse(1,Ny-3)];
            
            sex      = -spdiags(ax.Tsex2(x,3:Ny-1)',0,Ny-3,Ny-3);%-Wp.site.u_Inf*dxx(1:Ny-3,1:Ny-3)/2;    %v_{I,j}
            nex      = spdiags(ax.Tnex2(x,2:Ny-2)',0,Ny-3,Ny-3);%+Wp.site.u_Inf*dxx(1:Ny-3,1:Ny-3)/2;    %v_{I,j+1}
            Bx_E     = [sparse(1,Ny-3); sex] + [nex; sparse(1,Ny-3)];
            
            Axo((x-3)*(Ny-2)+1:(x-2)*(Ny-2) ,(x-3)*(Ny-3)+1:(x-1)*(Ny-3)) = [Bx_W Bx_E];
        end
        output.Axo = Axo;
        
        % For u-momentum equation
        ax.Tnx3D              = zeros(Nx,Ny);
        ax.Tsx3D              = zeros(Nx,Ny);
        
        ax.Tnx3D(1:Nx,1:Ny-1) = Rho/(Wp.turbine.Drotor)*dyy2(1:Nx,1:Ny-1).*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,1:Ny-1).^2))...
            .*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1));
        ax.Tsx3D(1:Nx,2:Ny)   = Rho/(Wp.turbine.Drotor)*dyy2(1:Nx,2:Ny).*(mixing_length(1:Nx,2:Ny).^2).*(dxx(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2))...
            .*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny));
       
        ax.aN             = ax.aN + ax.Tnx3D;
        ax.aS             = ax.aS + ax.Tsx3D;
        ax.aP             = ax.aP + ax.Tnx3D + ax.Tsx3D;
           
        % For v-momentum equation
        ay.Tey   = zeros(Nx,Ny);
        ay.Twy   = zeros(Nx,Ny);
        ay.Tny   = zeros(Nx,Ny);
        ay.Tsy   = zeros(Nx,Ny);
        
        ay.Tney  = zeros(Nx,Ny);
        ay.Tnwy  = zeros(Nx,Ny);
        
        %3
        ay.Tey(2:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx-1,1:Ny-1).^2).*(1./(dxx(1:Nx-1,1:Ny-1)))...
            .*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
        ay.Twy(1:Nx-1,2:Ny) = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(1./(dxx(1:Nx-1,1:Ny-1)))...
            .*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1));
        
        ay.aE             = ay.aE + ay.Tey;
        ay.aW             = ay.aW + ay.Twy;
        ay.aP             = ay.aP + ay.Tey + ay.Twy;
        
        ay.Tney(1:Nx-1,2:Ny) = 0.5*Rho*(mixing_length(1:Nx-1,1:Ny-1).^2).*(1./(dyy(2:Nx,1:Ny-1)))...
            .*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
        ay.Tnwy(2:Nx,2:Ny) = 0.5*Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(1./(dyy(2:Nx,1:Ny-1)))...
            .*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1));
        
        Ayo    = sparse((Ny-3)*(Nx-2),(Ny-2)*(Nx-3));
        By_W   = sparse(Ny-3,Ny-2); By_E   = sparse(Ny-3,Ny-2);
        
        for y= 3:Nx-1;
            swy      = spdiags(ay.Tnwy(y,3:Ny-1)',0,Ny-3,Ny-3);      %u_{i,J-1}
            nwy      = -spdiags(ay.Tnwy(y,3:Ny-1)',0,Ny-3,Ny-3);     %u_{i,J}
            By_W     = [swy sparse(Ny-3,1)] + [sparse(Ny-3,1) nwy];
            
            sey      = -spdiags(ay.Tney(y-1,3:Ny-1)',0,Ny-3,Ny-3);   %u_{i+1,J-1}
            ney      = spdiags(ay.Tney(y-1,3:Ny-1)',0,Ny-3,Ny-3);    %u_{i+1,J}
            By_E     = [sey sparse(Ny-3,1)] + [sparse(Ny-3,1) ney];
            
            Ayo((y-3)*(Ny-3)+1:(y-1)*(Ny-3) ,(y-3)*(Ny-2)+1:(y-2)*(Ny-2)) = [By_E; By_W];
        end
        output.Ayo = Ayo;
        
        
        %4
        ay.Tny(2:Nx,1:Ny-1) = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(dxx2(1:Nx-1,1:Ny-1)./(dyy(1:Nx-1,1:Ny-1).*dyy2(1:Nx-1,1:Ny-1)))...
            .*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
        ay.Tsy(1:Nx-1,2:Ny) = Rho*(mixing_length(2:Nx,2:Ny).^2).*(dxx2(1:Nx-1,1:Ny-1)./(dyy(1:Nx-1,1:Ny-1).*dyy2(1:Nx-1,1:Ny-1)))...
            .*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1));
        
        
        ay.aN             = ay.aN + ay.Tny;
        ay.aS             = ay.aS + ay.Tsy;
        ay.aP             = ay.aP + ay.Tsy + ay.Tny;
        
    case lower('WFSim4')
        
        % For u-momentum equation
        ax.Tex              = zeros(Nx,Ny);
        ax.Twx              = zeros(Nx,Ny);
        ax.Tnx              = zeros(Nx,Ny);
        ax.Tsx              = zeros(Nx,Ny);
        
        ax.Tnex             = zeros(Nx,Ny);
        ax.Tsex             = zeros(Nx,Ny);
        
        %1
        ax.Tex(1:Nx-1,2:Ny-1) = Rho*(mixing_length(1:Nx-1,2:Ny-1).^2).*(dyy2(1:Nx-1,2:Ny-1)./(dyy(1:Nx-1,2:Ny-1).*dxx2(1:Nx-1,2:Ny-1)))...
            .*abs(u(1:Nx-1,2:Ny-1)-u(1:Nx-1,1:Ny-2));
        ax.Twx(2:Nx,2:Ny-1)   = Rho*(mixing_length(2:Nx,2:Ny-1).^2).*(dyy2(1:Nx-1,2:Ny-1)./(dyy(2:Nx,2:Ny-1).*dxx2(2:Nx,2:Ny-1)))...
            .*abs(u(1:Nx-1,1:Ny-2)-u(1:Nx-1,2:Ny-1));
        
        ax.aE             = ax.aE + ax.Tex;
        ax.aW             = ax.aW + ax.Twx;
        ax.aP             = ax.aP + ax.Twx + ax.Tex;
        
        %2
        ax.Tnx(1:Nx,1:Ny-1)= Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,1:Ny-1).^2))...
            .*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1));
        ax.Tsx(1:Nx,2:Ny)  = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2))...
            .*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny));
        
        ax.aN             = ax.aN + ax.Tnx;
        ax.aS             = ax.aS + ax.Tsx;
        ax.aP             = ax.aP + ax.Tnx + ax.Tsx;
        
        ax.Tnex(1:Nx,1:Ny-1)= .5*Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(1./(dyy(1:Nx,1:Ny-1)))...
            .*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1));
        ax.Tsex(1:Nx,2:Ny)  = .5*Rho*(mixing_length(1:Nx,2:Ny).^2).*(1./(dyy(1:Nx,2:Ny)))...
            .*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny));
        
        Axo  = sparse((Ny-2)*(Nx-3),(Nx-2)*(Ny-3));
        Bx_W = sparse(Ny-2,Ny-3); Bx_E   = sparse(Ny-2,Ny-3);
        
        for x= 3:Nx-1;
            swx      = spdiags(ax.Tsex(x,3:Ny-1)',0,Ny-3,Ny-3);    %v_{I-1,j}
            nwx      = -spdiags(ax.Tnex(x,2:Ny-2)',0,Ny-3,Ny-3);    %v_{I-1,j+1}
            
            Bx_W     = [sparse(1,Ny-3); swx] + [nwx; sparse(1,Ny-3)];
            
            sex      = -spdiags(ax.Tsex(x,3:Ny-1)',0,Ny-3,Ny-3);    %v_{I,j}
            nex      = spdiags(ax.Tnex(x,2:Ny-2)',0,Ny-3,Ny-3);    %v_{I,j+1}
            Bx_E     = [sparse(1,Ny-3); sex] + [nex; sparse(1,Ny-3)];
            
            Axo((x-3)*(Ny-2)+1:(x-2)*(Ny-2) ,(x-3)*(Ny-3)+1:(x-1)*(Ny-3)) = [Bx_W Bx_E];
        end
        output.Axo = Axo;
        
        
        % For v-momentum equation
        ay.Tey   = zeros(Nx,Ny);
        ay.Twy   = zeros(Nx,Ny);
        ay.Tny   = zeros(Nx,Ny);
        ay.Tsy   = zeros(Nx,Ny);
        
        ay.Tney  = zeros(Nx,Ny);
        ay.Tnwy  = zeros(Nx,Ny);
        
        %3
        ay.Tey(2:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx-1,1:Ny-1).^2).*(1./(dxx(1:Nx-1,1:Ny-1)))...
            .*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
        ay.Twy(1:Nx-1,2:Ny) = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(1./(dxx(1:Nx-1,1:Ny-1)))...
            .*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1));
        
        ay.aE             = ay.aE + ay.Tey;
        ay.aW             = ay.aW + ay.Twy;
        ay.aP             = ay.aP + ay.Tey + ay.Twy;
        
        ay.Tney(1:Nx-1,2:Ny) = 0.5*Rho*(mixing_length(1:Nx-1,1:Ny-1).^2).*(1./(dyy(2:Nx,1:Ny-1)))...
            .*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
        ay.Tnwy(2:Nx,2:Ny) = 0.5*Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(1./(dyy(2:Nx,1:Ny-1)))...
            .*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1));
        
        Ayo    = sparse((Ny-3)*(Nx-2),(Ny-2)*(Nx-3));
        By_W   = sparse(Ny-3,Ny-2); By_E   = sparse(Ny-3,Ny-2);
        
        for y= 3:Nx-1;
            swy      = spdiags(ay.Tnwy(y,3:Ny-1)',0,Ny-3,Ny-3);      %u_{i,J-1}
            nwy      = -spdiags(ay.Tnwy(y,3:Ny-1)',0,Ny-3,Ny-3);     %u_{i,J}
            By_W     = [swy sparse(Ny-3,1)] + [sparse(Ny-3,1) nwy];
            
            sey      = -spdiags(ay.Tney(y-1,3:Ny-1)',0,Ny-3,Ny-3);   %u_{i+1,J-1}
            ney      = spdiags(ay.Tney(y-1,3:Ny-1)',0,Ny-3,Ny-3);    %u_{i+1,J}
            By_E     = [sey sparse(Ny-3,1)] + [sparse(Ny-3,1) ney];
            
            Ayo((y-3)*(Ny-3)+1:(y-1)*(Ny-3) ,(y-3)*(Ny-2)+1:(y-2)*(Ny-2)) = [By_E; By_W];
        end
        output.Ayo = Ayo;
        
        
        %4
        ay.Tny(2:Nx,1:Ny-1) = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(dxx2(1:Nx-1,1:Ny-1)./(dyy(1:Nx-1,1:Ny-1).*dyy2(1:Nx-1,1:Ny-1)))...
            .*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
        ay.Tsy(1:Nx-1,2:Ny) = Rho*(mixing_length(2:Nx,2:Ny).^2).*(dxx2(1:Nx-1,1:Ny-1)./(dyy(1:Nx-1,1:Ny-1).*dyy2(1:Nx-1,1:Ny-1)))...
            .*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1));
        
        
        ay.aN             = ay.aN + ay.Tny;
        ay.aS             = ay.aS + ay.Tsy;
        ay.aP             = ay.aP + ay.Tsy + ay.Tny;
        
    case lower('WFSim3')
        
        % For u-momentum equation
        ax.Tnx              = zeros(Nx,Ny);
        ax.Tsx              = zeros(Nx,Ny);
        
        ax.Tnx(2:Nx,1:Ny-1) = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(dxx(2:Nx,1:Ny-1)./(dyy(2:Nx,2:Ny).^2)).*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
        ax.Tsx(1:Nx-1,2:Ny) = Rho*(mixing_length(1:Nx-1,2:Ny).^2).*(dxx(2:Nx,2:Ny)./(dyy(2:Nx,2:Ny).^2)).*abs(u(2:Nx,1:Ny-1)-u(2:Nx,2:Ny));
        
        ax.aN             = ax.aN + ax.Tnx;
        ax.aS             = ax.aS + ax.Tsx;
        ax.aP             = ax.aP + ax.Tnx + ax.Tsx;
        
        % For v-momentum equation
        ay.Tey            = zeros(Nx,Ny);
        ay.Twy            = zeros(Nx,Ny);
        
        ay.Tey(1:Nx-1,1:Ny) = Rho*(mixing_length(1:Nx-1,1:Ny).^2).*(dyy(1:Nx-1,1:Ny)./(dxx(1:Nx-1,1:Ny).^2)).*abs(v(2:Nx,1:Ny)-v(1:Nx-1,1:Ny));
        ay.Twy(2:Nx,1:Ny)   = Rho*(mixing_length(2:Nx,1:Ny).^2).*(dyy(2:Nx,1:Ny)./(dxx(2:Nx,1:Ny).^2)).*abs(v(1:Nx-1,1:Ny)-v(2:Nx,1:Ny));
        
        ay.aE             = ay.aE + ay.Tey;
        ay.aW             = ay.aW + ay.Twy;
        ay.aP             = ay.aP + ay.Tey + ay.Twy;
        
        
    case lower('WFSim2')
        
        % For u-momentum equation
        ax.Tnx           = zeros(Nx,Ny);
        ax.Tsx           = zeros(Nx,Ny);
        
        ax.Tnx(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1));
        ax.Tsx(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny));
        
        ax.aN             = ax.aN + ax.Tnx;
        ax.aS             = ax.aS + ax.Tsx;
        ax.aP             = ax.aP + ax.Tnx + ax.Tsx;
        
        % For v-momentum equation
        ay.Tey = zeros(Nx,Ny);
        ay.Twy = zeros(Nx,Ny);
        
        ay.Tey(1:Nx-1,2:Ny)   = Rho*(mixing_length(1:Nx-1,2:Ny).^2)./dyy(1:Nx-1,2:Ny).*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
        ay.Twy(2:Nx,2:Ny)     = Rho*(mixing_length(2:Nx,2:Ny).^2)./dyy(1:Nx-1,2:Ny).*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1));
        
        % Define Ayo
        Ayo    = sparse((Ny-3)*(Nx-2),(Ny-2)*(Nx-3));
        By_W   = sparse(Ny-3,Ny-2); By_E   = sparse(Ny-3,Ny-2);
        
        for y= 3:Nx-1;
            swy      = spdiags(ay.Twy(y,3:Ny-1)',0,Ny-3,Ny-3); %u_{i,J-1}
            nwy      = -spdiags(ay.Twy(y,3:Ny-1)',0,Ny-3,Ny-3); %u_{i,J}
            By_W     = [swy sparse(Ny-3,1)] + [sparse(Ny-3,1) nwy];
            
            sey      = -spdiags(ay.Tey(y-1,3:Ny-1)',0,Ny-3,Ny-3);%u_{i+1,J-1}
            ney      = spdiags(ay.Tey(y-1,3:Ny-1)',0,Ny-3,Ny-3);%u_{i+1,J}
            By_E     = [sey sparse(Ny-3,1)] + [sparse(Ny-3,1) ney];
            
            Ayo((y-3)*(Ny-3)+1:(y-1)*(Ny-3) ,(y-3)*(Ny-2)+1:(y-2)*(Ny-2)) = [By_E; By_W];
        end
        
        output.Ayo = Ayo;
        
    case lower('WFSim1')
        
        % For u-momentum equation
        ax.Tnx              = zeros(Nx,Ny);
        ax.Tsx              = zeros(Nx,Ny);
        
        ax.Tnx(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2))...
            .*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1));
        ax.Tsx(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2))...
            .*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny));
                
        ax.aN             = ax.aN + ax.Tnx;
        ax.aS             = ax.aS + ax.Tsx;
        ax.aP             = ax.aP + ax.Tnx + ax.Tsx;
        
        % For v-momentum equation
        ay.Tey            = zeros(Nx,Ny);
        ay.Twy            = zeros(Nx,Ny);
        
        ay.Tey(1:Nx-1,1:Ny) = Rho*(mixing_length(1:Nx-1,1:Ny).^2).*(dyy(1:Nx-1,1:Ny)./(dxx(1:Nx-1,1:Ny).^2))...
            .*abs(v(2:Nx,1:Ny)-v(1:Nx-1,1:Ny));
        ay.Twy(2:Nx,1:Ny)   = Rho*(mixing_length(2:Nx,1:Ny).^2).*(dyy(2:Nx,1:Ny)./(dxx(2:Nx,1:Ny).^2))...
            .*abs(v(1:Nx-1,1:Ny)-v(2:Nx,1:Ny));
        
        ay.aE             = ay.aE + ay.Tey;
        ay.aW             = ay.aW + ay.Twy;
        ay.aP             = ay.aP + ay.Tey + ay.Twy;
        
        if Linearversion
            
            % For u-momentum equation
            %dTsxd1           = zeros(Nx,Ny);
            %dTsxd2           = zeros(Nx,Ny);
            %dTnxd1           = zeros(Nx,Ny);
            %dTnxd2           = zeros(Nx,Ny);
            
            % dTsx/du_(i,J)
            %dTsxd1(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx2(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*sign((u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny)));
            % dTsx/du_(i,J-1)
            %dTsxd2(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx2(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*-sign((u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny)));
            % dTnx/du_(i,J)
            %dTnxd1(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx2(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*-sign((u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1)));
            % dTnx/du_(i,J+1)
            %dTnxd2(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx2(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*sign((u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1)));
            
            dax.S(1:Nx,2:Ny)   = dax.S(1:Nx,2:Ny)   + ax.Tsx(1:Nx,2:Ny) ;
            dax.N(1:Nx,1:Ny-1) = dax.N(1:Nx,1:Ny-1) + ax.Tnx(1:Nx,1:Ny-1);
            dax.P(1:Nx,1:Ny-1) = dax.P(1:Nx,1:Ny-1) + ax.Tnx(1:Nx,1:Ny-1) + ax.Tsx(1:Nx,1:Ny-1) ;
            
            %dax.S(1:Nx,2:Ny)   = dax.S(1:Nx,2:Ny)+dTsxd2(1:Nx,2:Ny).*u(1:Nx,2:Ny)-dTsxd2(1:Nx,2:Ny).*u(1:Nx,1:Ny-1);
            %dax.N(1:Nx,1:Ny-1) = dax.N(1:Nx,1:Ny-1)-dTnxd2(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1)+dTnxd2(1:Nx,1:Ny-1).*u(1:Nx,2:Ny);
            %dax.P(1:Nx,1:Ny-1) = dax.P(1:Nx,1:Ny-1)+dTnxd1(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1)+dTsxd1(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1) ...
            %-dTnxd1(1:Nx,1:Ny-1).*u(1:Nx,2:Ny)-dTsxd1(1:Nx,1:Ny-1).*u(1:Nx,2:Ny);
            
            % For v-momentum equation
            day.E(1:Nx,2:Ny)   = day.E(1:Nx,2:Ny)   + ay.Tey(1:Nx,2:Ny) ;
            day.W(1:Nx,1:Ny-1) = day.W(1:Nx,1:Ny-1) + ay.Twy(1:Nx,1:Ny-1);
            day.P(1:Nx,1:Ny-1) = day.P(1:Nx,1:Ny-1) + ay.Tey(1:Nx,1:Ny-1) + ay.Twy(1:Nx,1:Ny-1) ;
        end;      
end
