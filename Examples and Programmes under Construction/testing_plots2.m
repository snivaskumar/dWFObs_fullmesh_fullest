clear all
close all
clc

%DExKF (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_ExKF_uinf5_noest_Nofus_turbC/workspace.mat')
% %DExKF (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_ExKF/workspace.mat')
% %DDExKF:2D (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_nofus/workspace.mat')
% %DDExKF:2D, turbC (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_nofus_turbC/workspace.mat')


load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_CIc0p2_turbC/workspace.mat')

for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u05{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE05{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError05{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow05{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp05{j} = d_Wp{j}; sol_array05{j} = d_sol_array{j}; sys05{j} = d_sys{j};
    scriptOptions05{j} = d_scriptOptions{j}; strucObs05{j} = d_strucObs{j};
end
plotdWFObs( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
% plotdCOV( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
plotdWFObs( Wp05{2},sol_array05{2},sys05{2},scriptOptions05{2},strucObs05{2} );

n = 2137;
P{1} = strucObs05{1}.Pk;
P{2} = strucObs05{2}.Pk;
x{1} = sol_array05{1}(end).x;
x{2} = sol_array05{2}(end).x;
d_x{1} = [1:2137];
d_x{2} = [1:2137];
filter = 'exkf';
[x,P] = fuse2(x{1},x{2},P{1},P{2},d_x{1}',d_x{2}','ICI','OPTIMAL',0.5,10,filter); 
% [x,P] = fuze(x,P,d_x,2,n,0,0,[1:n]','C',1,'CONSTANT'); % OPTIMAL % CONSTANT
strucObs05{1}.Pk = P;
strucObs05{2}.Pk = P;
sol_array05{1}(end).x = x;
sol_array05{2}(end).x = x;
Nx = Wp05{1}.mesh.Nx;
Ny = Wp05{1}.mesh.Ny;
u = 5*ones(Nx,Ny);
u(3:end-1,2:end-1)= reshape(x(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
[u,~,~] = Updateboundaries(Nx,Ny,u,u,u);
sol_array05{1}(end).u = u;
sol_array05{2}(end).u = u;
plotdWFObs( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
plotdCOV( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
% plotdWFObs( Wp05{2},sol_array05{2},sys05{2},scriptOptions05{2},strucObs05{2} );