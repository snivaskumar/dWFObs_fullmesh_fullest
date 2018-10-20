% function [Ix,Iy] = address();


load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4/workspace.mat')
Wp22 = Wp; sol_array22 = sol_array; sys22 = sys; 
scriptOptions22 = scriptOptions; strucObs22 = strucObs;
Nx  = Wp22.mesh.Nx;
Ny  = Wp22.mesh.Ny;
P = diag(strucObs22.Pk);
tmp                     = max(max(P));
Puv                     = tmp*ones(Nx,Ny);
Puv(3:end-1,2:end-1)    = reshape(P(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
[Puv,~,~] = Updateboundaries(Nx,Ny,Puv,Puv,Puv);    
[Ix,Iy] = find(Puv==5);

clearvars -except Ix Iy ex_Wp ex_sol_array ex_sys ex_scriptOptions ex_strucObs

% for i = 1:length(Ix)
%     Puv(Ix(i),Iy(i)) = 5;
% end