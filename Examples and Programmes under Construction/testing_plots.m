%% U_Inf = 11m/s and it is not estimated (ExKF: No Fusion, Fusion)

clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_sim_uinf11_noest/workspace.mat')
for i = 1:Wp.sim.NN
    u0(i)           = sol_array(i).site.u_Inf;
    RMSE0(i)        = sol_array(i).score.RMSE_cline;
    maxError0(i)    = sol_array(i).score.maxError;
    RMSE_flow0(i)   = sol_array(i).score.RMSE_flow;
    time0(i)        = sol_array(i).score.CPUtime;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys;
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NLInf/workspace.mat')
for i = 1:Wp.sim.NN
    u1(i)           = sol_array(i).site.u_Inf;
    RMSE1(i)        = sol_array(i).score.RMSE_cline;
    maxError1(i)    = sol_array(i).score.maxError;
    RMSE_flow1(i)   = sol_array(i).score.RMSE_flow;
    time1(i)        = sol_array(i).score.CPUtime;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys;
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%EnKF: 1D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_EnKF_1D_uinf11_noest_Qu0p1Qv0p1R0p1_en450_HS/workspace.mat')
for i = 1:Wp.sim.NN
    u2(i)           = sol_array(i).site.u_Inf;
    RMSE2(i)        = sol_array(i).score.RMSE_cline;
    maxError2(i)    = sol_array(i).score.maxError;
    RMSE_flow2(i)   = sol_array(i).score.RMSE_flow;
    time2(i)        = sol_array(i).score.CPUtime;
end
Wp2 = Wp; sol_array2 = sol_array; sys2 = sys;
scriptOptions2 = scriptOptions; strucObs2 = strucObs;
%DExKF: 2D (Arch.: 1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_2D_uinf11_noest_ICI0p5_Qu0p1Qv0p1R0p1_P20/workspace.mat')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_2D_uinf11_noest_no_Qu0p1Qv0p1R0p1_P20/workspace.mat')
for i = 1:Wp.sim.NN
    u4(i)           = sol_array(i).site.u_Inf;
    RMSE4(i)        = sol_array(i).score.RMSE_cline;
    maxError4(i)    = sol_array(i).score.maxError;
    RMSE_flow4(i)   = sol_array(i).score.RMSE_flow;
    time4(i)        = sol_array(i).score.CPUtime;
end
Wp4 = Wp; sol_array4 = sol_array; sys4 = sys;
scriptOptions4 = scriptOptions; strucObs4 = strucObs;
%ExKF: No Fusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/turb2axialm_ExKF_uinf11_noest_nofus_turbC_Qu0p1Qv0p1R0p1/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u01{j}(i)           = d_sol_array{j}(i).site.u_Inf;
        RMSE01{j}(i)        = d_sol_array{j}(i).score.RMSE_cline;
        maxError01{j}(i)    = d_sol_array{j}(i).score.maxError;
        RMSE_flow01{j}(i)   = d_sol_array{j}(i).score.RMSE_flow;
        time01{j}(i)        = d_sol_array{j}(i).score.CPUtime;
    end
    tim01(j) = sum(time01{j})/d_Wp{j}.sim.NN;
    Wp01{j} = d_Wp{j}; sol_array01{j} = d_sol_array{j}; sys01{j} = d_sys{j};
    scriptOptions01{j} = d_scriptOptions{j}; strucObs01{j} = d_strucObs{j};
end
%ExKF: Fusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/turb2axialm_ExKF_uinf11_noest_ICI0p5_turbC_Qu0p1Qv0p1R0p1/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u02{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE02{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError02{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow02{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
        time02{j}(i)        = d_sol_array{j}(i).score.CPUtime;
    end
    tim02(j) = sum(time02{j})/d_Wp{j}.sim.NN;
    Wp02{j} = d_Wp{j}; sol_array02{j} = d_sol_array{j}; sys02{j} = d_sys{j};
    scriptOptions02{j} = d_scriptOptions{j}; strucObs02{j} = d_strucObs{j};
end
%DExKF: No Fusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/turb2axialm_DExKF_2D_uinf11_noest_nofus_turbC_Qu0p1Qv0p1R0p1/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u03{j}(i)           = d_sol_array{j}(i).site.u_Inf;
        RMSE03{j}(i)        = d_sol_array{j}(i).score.RMSE_cline;
        maxError03{j}(i)    = d_sol_array{j}(i).score.maxError;
        RMSE_flow03{j}(i)   = d_sol_array{j}(i).score.RMSE_flow;
        time03{j}(i)        = d_sol_array{j}(i).score.CPUtime;
    end
    tim03(j) = sum(time03{j})/d_Wp{j}.sim.NN;
    Wp03{j} = d_Wp{j}; sol_array03{j} = d_sol_array{j}; sys03{j} = d_sys{j};
    scriptOptions03{j} = d_scriptOptions{j}; strucObs03{j} = d_strucObs{j};
end
%DExKF: Fusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/turb2axialm_DExKF_2D_uinf11_noest_ICI0p5_turbC_Qu0p1Qv0p1R0p1/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u04{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE04{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError04{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow04{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
        time04{j}(i)        = d_sol_array{j}(i).score.CPUtime;
    end
    tim04(j) = sum(time04{j})/d_Wp{j}.sim.NN;
    Wp04{j} = d_Wp{j}; sol_array04{j} = d_sol_array{j}; sys04{j} = d_sys{j};
    scriptOptions04{j} = d_scriptOptions{j}; strucObs04{j} = d_strucObs{j};
end


%% U_Inf = 11m/s and it is not estimated (SIM, EnKF, DDExKF_IFAC:1D,2D)

clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_sim_uinf11_noest/workspace.mat')
for i = 1:Wp.sim.NN
    time0(i) = sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys;
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_uinf11_noest_Qev0p5/workspace.mat')
for i = 1:Wp.sim.NN
time1(i) = sol_array(i).score.CPUtime;
RMSE1(i) = sol_array(i).score.RMSE_cline;
maxError1(i) = sol_array(i).score.maxError;
RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys;
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%EnKF: 2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_2D_uinf11_noest_Qev0p5/workspace.mat')
for i = 1:Wp.sim.NN
time1_2D(i) = sol_array(i).score.CPUtime;
RMSE1_2D(i) = sol_array(i).score.RMSE_cline;
maxError1_2D(i) = sol_array(i).score.maxError;
RMSE_flow1_2D(i) = sol_array(i).score.RMSE_flow;
end
Wp1_2D = Wp; sol_array1_2D = sol_array; sys1_2D = sys;
scriptOptions1_2D = scriptOptions; strucObs1_2D = strucObs;
%EnKF: 3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_3D_uinf11_noest_Qev0p5/workspace.mat')
for i = 1:Wp.sim.NN
time1_3D(i) = sol_array(i).score.CPUtime;
RMSE1_3D(i) = sol_array(i).score.RMSE_cline;
maxError1_3D(i) = sol_array(i).score.maxError;
RMSE_flow1_3D(i) = sol_array(i).score.RMSE_flow;
end
Wp1_3D = Wp; sol_array1_3D = sol_array; sys1_3D = sys;
scriptOptions1_3D = scriptOptions; strucObs1_3D = strucObs;
%EnKF: 4D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_4D_uinf11_noest_Qev0p5/workspace.mat')
for i = 1:Wp.sim.NN
time1_4D(i) = sol_array(i).score.CPUtime;
RMSE1_4D(i) = sol_array(i).score.RMSE_cline;
maxError1_4D(i) = sol_array(i).score.maxError;
RMSE_flow1_4D(i) = sol_array(i).score.RMSE_flow;
end
Wp1_4D = Wp; sol_array1_4D = sol_array; sys1_4D = sys;
scriptOptions1_4D = scriptOptions; strucObs1_4D = strucObs;

figure, plot(RMSE_flow1), hold on, plot(RMSE_flow1_2D), plot(RMSE_flow1_3D), plot(RMSE_flow1_4D)
legend('EnKF: 1D','EnKF: 2D','EnKF: 3D','EnKF: 4D')
title('Flow Error')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.15);

figure, plot(RMSE1), hold on, plot(RMSE1_2D), plot(RMSE1_3D), plot(RMSE1_4D)
legend('EnKF: 1D','EnKF: 2D','EnKF: 3D','EnKF: 4D')
title('Centreline Error')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.15);
%ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest/workspace.mat')
for i = 1:Wp.sim.NN
    time2(i) = sol_array(i).score.CPUtime;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
    maxError2(i) = sol_array(i).score.maxError;
    RMSE_flow2(i) = sol_array(i).score.RMSE_flow;
end
Wp2 = Wp; sol_array2 = sol_array; sys2 = sys;
scriptOptions2 = scriptOptions; strucObs2 = strucObs;
%ExKF: 1D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_1D_uinf11_noest/workspace.mat')
for i = 1:Wp.sim.NN
    time2_1D(i) = sol_array(i).score.CPUtime;
    RMSE2_1D(i) = sol_array(i).score.RMSE_cline;
    maxError2_1D(i) = sol_array(i).score.maxError;
    RMSE_flow2_1D(i) = sol_array(i).score.RMSE_flow;
end
Wp2_1D = Wp; sol_array2_1D = sol_array; sys2_1D = sys;
scriptOptions2_1D = scriptOptions; strucObs2_1D = strucObs;
%ExKF: 2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_2D_uinf11_noest/workspace.mat')
for i = 1:Wp.sim.NN
    time2_2D(i) = sol_array(i).score.CPUtime;
    RMSE2_2D(i) = sol_array(i).score.RMSE_cline;
    maxError2_2D(i) = sol_array(i).score.maxError;
    RMSE_flow2_2D(i) = sol_array(i).score.RMSE_flow;
end
Wp2_2D = Wp; sol_array2_2D = sol_array; sys2_2D = sys;
scriptOptions2_2D = scriptOptions; strucObs2_2D = strucObs;
%ExKF: 3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_3D_uinf11_noest/workspace.mat')
for i = 1:Wp.sim.NN
time2_3D(i) = sol_array(i).score.CPUtime;
RMSE2_3D(i) = sol_array(i).score.RMSE_cline;
maxError2_3D(i) = sol_array(i).score.maxError;
RMSE_flow2_3D(i) = sol_array(i).score.RMSE_flow;
end
Wp2_3D = Wp; sol_array2_3D = sol_array; sys2_3D = sys;
scriptOptions2_3D = scriptOptions; strucObs2_3D = strucObs;
%ExKF: 4D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_4D_uinf11_noest/workspace.mat')
for i = 1:Wp.sim.NN
time2_4D(i) = sol_array(i).score.CPUtime;
RMSE2_4D(i) = sol_array(i).score.RMSE_cline;
maxError2_4D(i) = sol_array(i).score.maxError;
RMSE_flow2_4D(i) = sol_array(i).score.RMSE_flow;
end
Wp2_4D = Wp; sol_array2_4D = sol_array; sys2_4D = sys;
scriptOptions2_4D = scriptOptions; strucObs2_4D = strucObs;

figure, plot(RMSE_flow2), hold on, plot(RMSE_flow2_2D), plot(RMSE_flow2_3D), plot(RMSE_flow2_4D)
legend('ExKF: 1D','ExKF: 2D','ExKF: 3D','ExKF: 4D')
title('Flow Error')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.15);

figure, plot(RMSE2), hold on, plot(RMSE2_2D), plot(RMSE2_3D), plot(RMSE2_4D)
legend('ExKF: 1D','ExKF: 2D','ExKF: 3D','ExKF: 4D')
title('Centreline Error')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.15);
%IFAC_2DZE_extremeopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf11_noest_eopti_1em4/workspace.mat')
for i = 1:Wp.sim.NN
    time22(i) = sol_array(i).score.CPUtime;
    RMSE22(i) = sol_array(i).score.RMSE_cline;
    maxError22(i) = sol_array(i).score.maxError;
    RMSE_flow22(i) = sol_array(i).score.RMSE_flow;
end
Wp22 = Wp; sol_array22 = sol_array; sys22 = sys; 
scriptOptions22 = scriptOptions; strucObs22 = strucObs;
%DExKF: turbC (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_ExKF_uinf11_noest_Nofus_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u01{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE01{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError01{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow01{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp01{j} = d_Wp{j}; sol_array01{j} = d_sol_array{j}; sys01{j} = d_sys{j};
    scriptOptions01{j} = d_scriptOptions{j}; strucObs01{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (No Fusion), P_unest = 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf11_noest_nofus_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u02{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE02{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError02{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow02{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp02{j} = d_Wp{j}; sol_array02{j} = d_sol_array{j}; sys02{j} = d_sys{j};
    scriptOptions02{j} = d_scriptOptions{j}; strucObs02{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (CIc0p2), P0 = 1m/s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf11_noest_CIc0p2_turbC_Po1/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u03{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE03{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError03{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow03{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp03{j} = d_Wp{j}; sol_array03{j} = d_sol_array{j}; sys03{j} = d_sys{j};
    scriptOptions03{j} = d_scriptOptions{j}; strucObs03{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (CIc0p2), P_unest = 0.5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf11_noest_CIc0p2_turbC_P0p5/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u04{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE04{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError04{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow04{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp04{j} = d_Wp{j}; sol_array04{j} = d_sol_array{j}; sys04{j} = d_sys{j};
    scriptOptions04{j} = d_scriptOptions{j}; strucObs04{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (CIc0p2), P_unest = 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf11_noest_CIc0p2_turbC/workspace.mat')
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
%DDExKF:2D, turbC (CIc0p2), P_unest = 1e1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf11_noest_CIc0p2_turbC_P1e1/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u06{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE06{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError06{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow06{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp06{j} = d_Wp{j}; sol_array06{j} = d_sol_array{j}; sys06{j} = d_sys{j};
    scriptOptions06{j} = d_scriptOptions{j}; strucObs06{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (CIc0p2), P_unest = 1e2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf11_noest_CIc0p2_turbC_P1e2/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u07{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE07{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError07{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow07{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp07{j} = d_Wp{j}; sol_array07{j} = d_sol_array{j}; sys07{j} = d_sys{j};
    scriptOptions07{j} = d_scriptOptions{j}; strucObs07{j} = d_strucObs{j};
end
%DDExKF:3D, turbC (CIc0p2), P_unest = 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_3D_IFAC_uinf11_noest_CIc0p2_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u08{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE08{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError08{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow08{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp08{j} = d_Wp{j}; sol_array08{j} = d_sol_array{j}; sys08{j} = d_sys{j};
    scriptOptions08{j} = d_scriptOptions{j}; strucObs08{j} = d_strucObs{j};
end
%DDExKF:3D, turbC (CIc0p2), P_unest = 1e1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_3D_IFAC_uinf11_noest_CIc0p2_turbC_P1e1/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u09{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE09{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError09{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow09{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp09{j} = d_Wp{j}; sol_array09{j} = d_sol_array{j}; sys09{j} = d_sys{j};
    scriptOptions09{j} = d_scriptOptions{j}; strucObs09{j} = d_strucObs{j};
end
%DDExKF:3D, turbC (CIc0p2), P_unest = 1e2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_3D_IFAC_uinf11_noest_CIc0p2_turbC_P1e2/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u10{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE10{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError10{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow10{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp10{j} = d_Wp{j}; sol_array10{j} = d_sol_array{j}; sys10{j} = d_sys{j};
    scriptOptions10{j} = d_scriptOptions{j}; strucObs10{j} = d_strucObs{j};
end
%DDExKF:4D, turbC (ICIc0p5), P_unest = 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_4D_IFAC_uinf11_noest_ICIc0p5_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u11{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE11{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError11{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow11{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp11{j} = d_Wp{j}; sol_array11{j} = d_sol_array{j}; sys11{j} = d_sys{j};
    scriptOptions11{j} = d_scriptOptions{j}; strucObs11{j} = d_strucObs{j};
end
% tp05 = RMSE_flow05;   tp06 = RMSE_flow06;
% tp07 = RMSE_flow07;   tp08 = RMSE_flow08;
% tp09 = RMSE_flow09;   tp10 = RMSE_flow10;
% tp11 = RMSE_flow11;
% for j = 1:2
%     for i = 1:length(RMSE05{1})
%         if mean(RMSE05{1}) >= 3
%             tp05{j}(i) = RMSE_flow05{j}(i-1);
%         end
%     end
% end
% for j = 1:2
%     for i = 1:length(RMSE06{1})
%         if mean(RMSE06{1}) >= 3
%             tp06{j}(i) = RMSE_flow06{j}(i-1);
%         end
%     end
% end
% for j = 1:2
%     for i = 1:length(RMSE07{1})
%         if mean(RMSE07{1}) >= 3
%             tp07{j}(i) = RMSE_flow07{j}(i-1);
%         end
%     end
% end
% for j = 1:2
%     for i = 1:length(RMSE08{1})
%         if mean(RMSE08{1}) >= 3
%             tp08{j}(i) = RMSE_flow08{j}(i-1);
%         end
%     end
% end
% for j = 1:2
%     for i = 1:length(RMSE09{1})
%         if mean(RMSE09{1}) >= 3
%             tp09{j}(i) = RMSE_flow09{j}(i-1);
%         end
%     end
% end
% for j = 1:2
%     for i = 1:length(RMSE10{1})
%         if mean(RMSE10{1}) >= 3
%             tp10{j}(i) = RMSE_flow10{j}(i-1);
%         end
%     end
% end
% for j = 1:2
%     for i = 1:length(RMSE11{1})
%         if mean(RMSE11{1}) >= 3
%             tp11{j}(i) = RMSE_flow11{j}(i-1);
%         end
%     end
% end

% figure, plot(RMSE0), hold on, plot(RMSE1), plot(RMSE22), 
% plot(RMSE03{1}),plot(RMSE04{1}),
plot(RMSE05{1}), hold on, plot(RMSE06{1}),plot(RMSE07{1}),
plot(RMSE08{1}),plot(RMSE09{1}),plot(RMSE10{1}),
plot(RMSE11{1})
title('Centreline Error')
% legend('sim','EnKF','DExKF: 2D (arch. 1)',...
%     'DDExKF:2D (CI_P0_1)','DDExKF:2D (CI_Pun_0.5)',...
%     'DDExKF:2D (CI_Pun_5)','DDExKF:2D (CI_Pun_1e1)','DDExKF:2D (CI_Pun_1e2)',...
%     'DDExKF:3D (CI_Pun_5)','DDExKF:3D (CI_Pun_1e1)','DDExKF:3D (CI_Pun_1e2)',...
%     'DDExKF:4D (CI_Pun_5)'...
%     );
legend('DDExKF:2D (CI_Pun_5)','DDExKF:2D (CI_Pun_1e1)','DDExKF:2D (CI_Pun_1e2)',...
    'DDExKF:3D (CI_Pun_5)','DDExKF:3D (CI_Pun_1e1)','DDExKF:3D (CI_Pun_1e2)',...
    'DDExKF:4D (CI_Pun_5)'...
    );
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.15);
% figure, plot(RMSE_flow0), hold on, plot(RMSE_flow1),plot(RMSE_flow22),
% plot(RMSE_flow03{1}),plot(RMSE_flow04{1}),
plot(RMSE_flow05{1}), hold on, plot(RMSE_flow06{1}),plot(RMSE_flow07{1}),
plot(RMSE_flow08{1}),plot(RMSE_flow09{1}),plot(RMSE_flow10{1}),
plot(RMSE_flow11{1})
title('Flow Error')
legend('DDExKF:2D (CI_Pun_5)','DDExKF:2D (CI_Pun_1e1)','DDExKF:2D (CI_Pun_1e2)',...
    'DDExKF:3D (CI_Pun_5)','DDExKF:3D (CI_Pun_1e1)','DDExKF:3D (CI_Pun_1e2)',...
    'DDExKF:4D (CI_Pun_5)');
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.15);

%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
%EnKF_2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp1_2D,sol_array1_2D,sys1_2D,scriptOptions1_2D,strucObs1_2D );
%EnKF_3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp1_3D,sol_array1_3D,sys1_3D,scriptOptions1_3D,strucObs1_3D );
%EnKF_4D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp1_4D,sol_array1_4D,sys1_4D,scriptOptions1_4D,strucObs1_4D );
%ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp2,sol_array2,sys2,scriptOptions2,strucObs2 );
%ExKF_1D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp2_1D,sol_array2_1D,sys2_1D,scriptOptions2_1D,strucObs2_1D );
%ExKF_2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp2_2D,sol_array2_2D,sys2_2D,scriptOptions2_2D,strucObs2_2D );
%ExKF_3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp2_3D,sol_array2_3D,sys2_3D,scriptOptions2_3D,strucObs2_3D );
%ExKF_4D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotWFObs( Wp2_4D,sol_array2_4D,sys2_4D,scriptOptions2_4D,strucObs2_4D );
%IFAC_2DZE_extremeopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
plotWFObs( Wp22,sol_array22,sys22,scriptOptions22,strucObs22 );
%DExKF:2D, turbC (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp01{1},sol_array01{1},sys01{1},scriptOptions01{1},strucObs01{1} );
plotdWFObs( Wp01{2},sol_array01{2},sys01{2},scriptOptions01{2},strucObs01{2} );
%DDExKF:2D, turbC (No Fusion), P_unest = 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp02{1},sol_array02{1},sys02{1},scriptOptions02{1},strucObs02{1} );
plotdWFObs( Wp02{2},sol_array02{2},sys02{2},scriptOptions02{2},strucObs02{2} );
%DDExKF:2D, turbC (CIc0p2), P0 = 1m/s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp03{1},sol_array03{1},sys03{1},scriptOptions03{1},strucObs03{1} );
%DDExKF:2D, turbC (CIc0p2), P_unest = 0.5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp04{1},sol_array04{1},sys04{1},scriptOptions04{1},strucObs04{1} );
%DDExKF:2D, turbC (CIc0p2), P_unest = 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
%DDExKF:2D, turbC (CIc0p2), P_unest = 1e1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp06{1},sol_array06{1},sys06{1},scriptOptions06{1},strucObs06{1} );
%DDExKF:2D, turbC (CIc0p2), P_unest = 1e2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp07{1},sol_array07{1},sys07{1},scriptOptions07{1},strucObs07{1} );
%DDExKF:3D, turbC (CIc0p2), P_unest = 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp08{1},sol_array08{1},sys08{1},scriptOptions08{1},strucObs08{1} );
%DDExKF:3D, turbC (CIc0p2), P_unest = 1e1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp09{1},sol_array09{1},sys09{1},scriptOptions09{1},strucObs09{1} );
%DDExKF:3D, turbC (CIc0p2), P_unest = 1e2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp10{1},sol_array10{1},sys10{1},scriptOptions10{1},strucObs10{1} );
%DDExKF:4D, turbC (ICIc0p5), P_unest = 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdWFObs( Wp11{1},sol_array11{1},sys11{1},scriptOptions11{1},strucObs11{1} );

%% U_Inf = 5m/s and it is not estimated (SIM, EnKF, DDExKF_IFAC:1D,2D)

clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_sim_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    time0(i) = sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys;
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time1(i) = sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys;
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%IFAC_2DZE_extremeopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim22,timeWFSim22,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time22(i) = sol_array(i).score.CPUtime;
    RMSE22(i) = sol_array(i).score.RMSE_cline;
    maxError22(i) = sol_array(i).score.maxError;
    RMSE_flow22(i) = sol_array(i).score.RMSE_flow;
end
Wp22 = Wp; sol_array22 = sol_array; sys22 = sys; 
scriptOptions22 = scriptOptions; strucObs22 = strucObs;
% %DExKF (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_ExKF/workspace.mat')
% for j = 1:2
%     for i = 1:d_Wp{j}.sim.NN
%         u00{j}(i) = d_sol_array{j}(i).site.u_Inf;
%         RMSE00{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
%         maxError00{j}(i) = d_sol_array{j}(i).score.maxError;
%         RMSE_flow00{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
%     end
%     Wp00{j} = d_Wp{j}; sol_array00{j} = d_sol_array{j}; sys00{j} = d_sys{j};
%     scriptOptions00{j} = d_scriptOptions{j}; strucObs00{j} = d_strucObs{j};
% end
% %DExKF (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_ExKF_uinf5_noest_Nofus_turbC/workspace.mat')
% for j = 1:2
%     for i = 1:d_Wp{j}.sim.NN
%         u000{j}(i) = d_sol_array{j}(i).site.u_Inf;
%         RMSE000{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
%         maxError000{j}(i) = d_sol_array{j}(i).score.maxError;
%         RMSE_flow000{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
%     end
%     Wp000{j} = d_Wp{j}; sol_array000{j} = d_sol_array{j}; sys000{j} = d_sys{j};
%     scriptOptions000{j} = d_scriptOptions{j}; strucObs000{j} = d_strucObs{j};
% end
% %DDExKF:1D (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_1D_IFAC_uinf5_noest/workspace.mat')
% for j = 1:2
%     for i = 1:d_Wp{j}.sim.NN
%         u01{j}(i) = d_sol_array{j}(i).site.u_Inf;
%         RMSE01{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
%         maxError01{j}(i) = d_sol_array{j}(i).score.maxError;
%         RMSE_flow01{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
%     end
%     Wp01{j} = d_Wp{j}; sol_array01{j} = d_sol_array{j}; sys01{j} = d_sys{j};
%     scriptOptions01{j} = d_scriptOptions{j}; strucObs01{j} = d_strucObs{j};
% end
%DDExKF:2D (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u02{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE02{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError02{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow02{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp02{j} = d_Wp{j}; sol_array02{j} = d_sol_array{j}; sys02{j} = d_sys{j};
    scriptOptions02{j} = d_scriptOptions{j}; strucObs02{j} = d_strucObs{j};
end
% %DDExKF:3D (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_3D_IFAC_uinf5_noest/workspace.mat')
% for j = 1:2
%     for i = 1:d_Wp{j}.sim.NN
%         u03{j}(i) = d_sol_array{j}(i).site.u_Inf;
%         RMSE03{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
%         maxError03{j}(i) = d_sol_array{j}(i).score.maxError;
%         RMSE_flow03{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
%     end
%     Wp03{j} = d_Wp{j}; sol_array03{j} = d_sol_array{j}; sys03{j} = d_sys{j};
%     scriptOptions03{j} = d_scriptOptions{j}; strucObs03{j} = d_strucObs{j};
% end
% %DDExKF:2D:SQUARE (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2Dsq_IFAC_uinf5_noest/workspace.mat')
% for j = 1:2
%     for i = 1:d_Wp{j}.sim.NN
%         u04{j}(i) = d_sol_array{j}(i).site.u_Inf;
%         RMSE04{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
%         maxError04{j}(i) = d_sol_array{j}(i).score.maxError;
%         RMSE_flow04{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
%     end
%     Wp04{j} = d_Wp{j}; sol_array04{j} = d_sol_array{j}; sys04{j} = d_sys{j};
%     scriptOptions04{j} = d_scriptOptions{j}; strucObs04{j} = d_strucObs{j};
% end
%DDExKF:2D (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_nofus/workspace.mat')
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
%DDExKF:2D (CIN)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_CIN/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u06{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE06{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError06{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow06{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp06{j} = d_Wp{j}; sol_array06{j} = d_sol_array{j}; sys06{j} = d_sys{j};
    scriptOptions06{j} = d_scriptOptions{j}; strucObs06{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_nofus_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u07{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE07{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError07{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow07{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp07{j} = d_Wp{j}; sol_array07{j} = d_sol_array{j}; sys07{j} = d_sys{j};
    scriptOptions07{j} = d_scriptOptions{j}; strucObs07{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_IFAC_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u08{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE08{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError08{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow08{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp08{j} = d_Wp{j}; sol_array08{j} = d_sol_array{j}; sys08{j} = d_sys{j};
    scriptOptions08{j} = d_scriptOptions{j}; strucObs08{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (CIN)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_CIN_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u09{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE09{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError09{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow09{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp09{j} = d_Wp{j}; sol_array09{j} = d_sol_array{j}; sys09{j} = d_sys{j};
    scriptOptions09{j} = d_scriptOptions{j}; strucObs09{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (ICI)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_ICIc0p5_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u10{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE10{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError10{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow10{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp10{j} = d_Wp{j}; sol_array10{j} = d_sol_array{j}; sys10{j} = d_sys{j};
    scriptOptions10{j} = d_scriptOptions{j}; strucObs10{j} = d_strucObs{j};
end
%DDExKF:2D, turbNEW (ICI)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_ICI/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u11{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE11{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError11{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow11{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp11{j} = d_Wp{j}; sol_array11{j} = d_sol_array{j}; sys11{j} = d_sys{j};
    scriptOptions11{j} = d_scriptOptions{j}; strucObs11{j} = d_strucObs{j};
end
%DDExKF:2D, turbC (CIc0p2)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear d_scriptOptions d_sol_array d_strucObs d_sys d_Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/axi_2turb_alm_turb_DExKF_2D_IFAC_uinf5_noest_CIc0p2_turbC/workspace.mat')
for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u12{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE12{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError12{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow12{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp12{j} = d_Wp{j}; sol_array12{j} = d_sol_array{j}; sys12{j} = d_sys{j};
    scriptOptions12{j} = d_scriptOptions{j}; strucObs12{j} = d_strucObs{j};
end

figure, plot(RMSE0), hold on, plot(RMSE1), plot(RMSE22), 
% plot(RMSE07{1}),plot(RMSE07{2}),
% 'DDExKF:2D,sys1 (No Fusion_turbC)','DDExKF:2D,sys2 (No Fusion_turbC)',...
plot(RMSE08{1}),plot(RMSE09{1}),
plot(RMSE10{1}),plot(RMSE12{1})
title('Centreline Error')
legend('sim','EnKF','DExKF: 2D (arch. 1)',...
    'DDExKF:2D (IFAC_turbC)','DDExKF:2D (CIN_turbC)',...
    'DDExKF:2D (ICI_turbC)','DDExKF:2D (CI_turbC)')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.15);
figure, plot(RMSE_flow0), hold on, plot(RMSE_flow1),plot(RMSE_flow22),
% plot(RMSE_flow07{1}),plot(RMSE_flow07{2}),
plot(RMSE_flow08{1}),plot(RMSE_flow09{1}),
plot(RMSE_flow10{1}),plot(RMSE_flow12{1})
title('Flow Error')
legend('sim','EnKF','DExKF: 2D (arch. 1)',...
    'DDExKF:2D (IFAC_turbC)','DDExKF:2D (CIN_turbC)',...
    'DDExKF:2D (ICI_turbC)','DDExKF:2D (CI_turbC)')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.15);

% figure, plot(RMSE0), hold on, plot(RMSE1), plot(RMSE22), 
% % plot(RMSE02{1}),plot(RMSE05{1}),plot(RMSE05{2}),plot(RMSE06{1}),
% % 'DDExKF:2D (IFAC)','DDExKF:2D,sys1 (No Fusion)','DDExKF:2D,sys2 (No Fusion)','DDExKF:2D (CIN)',...
% % plot(RMSE07{1}),plot(RMSE07{2}),plot(RMSE08{1}),plot(RMSE09{1}),
% % 'DDExKF:2D,sys1 (No Fusion_turbC)','DDExKF:2D,sys2 (No Fusion_turbC)','DDExKF:2D (IFAC_turbC)','DDExKF:2D (CIN_turbC)',...
% plot(RMSE10{1}),plot(RMSE11{1}),plot(RMSE12{1})
% title('Centreline Error')
% legend('sim','EnKF','DExKF: 2D (arch. 1)',...
%     'DDExKF:2D (ICI_turbC)','DDExKF:2D (ICI_turbNEW)','DDExKF:2D (CI_turbC)')
% lq = findobj(gcf,'type','line');
% set(lq,'linewidth',1.15);
% figure, plot(RMSE_flow0), hold on, plot(RMSE_flow1),plot(RMSE_flow22),
% % plot(RMSE_flow02{1}),plot(RMSE_flow05{1}),plot(RMSE_flow05{2}),plot(RMSE_flow06{1}),
% % plot(RMSE_flow07{1}),plot(RMSE_flow07{2}),plot(RMSE_flow08{1}),plot(RMSE_flow09{1}),
% plot(RMSE_flow10{1}),plot(RMSE_flow11{1}),plot(RMSE_flow12{1})
% title('Flow Error')
% legend('sim','EnKF','DExKF: 2D (arch. 1)',...
%     'DDExKF:2D (ICI_turbC)','DDExKF:2D (ICI_turbNEW)','DDExKF:2D (CI_turbC)')
% lq = findobj(gcf,'type','line');
% set(lq,'linewidth',1.15);

% % %SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
% % EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
% %DExKF: 2D (archiecture 1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotWFObs( Wp22,sol_array22,sys22,scriptOptions22,strucObs22 );
% % %DExKF (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plotdWFObs( Wp00{1},sol_array00{1},sys00{1},scriptOptions00{1},strucObs00{1});
% % %DDExKF: 1D (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plotdWFObs( Wp01{1},sol_array01{1},sys01{1},scriptOptions01{1},strucObs01{1} );
% % DDExKF: 2D (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotdWFObs( Wp02{1},sol_array02{1},sys02{1},scriptOptions02{1},strucObs02{1} );
% % %DDExKF: 3D (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plotdWFObs( Wp03{1},sol_array03{1},sys03{1},scriptOptions03{1},strucObs03{1} );
% % %DDExKF: 2Dsq (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plotdWFObs( Wp04{1},sol_array04{1},sys04{1},scriptOptions04{1},strucObs04{1} );
% %DDExKF: 2D, turbNEW (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotdWFObs( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
% plotdWFObs( Wp05{2},sol_array05{2},sys05{2},scriptOptions05{2},strucObs05{2} );
% %DDExKF: 2D, turbNEW (CIN)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotdWFObs( Wp06{1},sol_array06{1},sys06{1},scriptOptions06{1},strucObs06{1} );
% %DDExKF:2D, turbC (No Fusion)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotdWFObs( Wp07{1},sol_array07{1},sys07{1},scriptOptions07{1},strucObs07{1} );
% plotdWFObs( Wp07{2},sol_array07{2},sys07{2},scriptOptions07{2},strucObs07{2} );
% %DDExKF:2D, turbC (IFAC)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotdWFObs( Wp08{1},sol_array08{1},sys08{1},scriptOptions08{1},strucObs08{1} );
% %DDExKF:2D, turbC (CIN)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotdWFObs( Wp09{1},sol_array09{1},sys09{1},scriptOptions09{1},strucObs09{1} );
% %DDExKF:2D, turbC (ICI)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotdWFObs( Wp10{1},sol_array10{1},sys10{1},scriptOptions10{1},strucObs10{1} );
% % %DDExKF:2D, turbNEW (ICI)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plotdWFObs( Wp11{1},sol_array11{1},sys11{1},scriptOptions11{1},strucObs11{1} );
% %DDExKF:2D, turbC (CIc0p2)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotdWFObs( Wp12{1},sol_array12{1},sys12{1},scriptOptions12{1},strucObs12{1} );
%% U_Inf = 11m/s and it is not estimated (EnKF, DExKF_IFAC:1D,2D)
% With E
clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_sim_uinf11_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_uinf_noest\workspace.mat')
for i = 1:Wp.sim.NN
    u0(i) = sol_array(i).site.u_Inf;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys; 
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF:1D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_uinf11_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_uinf_noest\workspace.mat')
for i = 1:Wp.sim.NN
    u1(i) = sol_array(i).site.u_Inf;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys; 
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%EnKF:2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_2D_uinf11_noest/workspace.mat')
for i = 1:Wp.sim.NN
    u2(i) = sol_array(i).site.u_Inf;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
    maxError2(i) = sol_array(i).score.maxError;
    RMSE_flow2(i) = sol_array(i).score.RMSE_flow;
end
Wp2 = Wp; sol_array2 = sol_array; sys2 = sys; 
scriptOptions2 = scriptOptions; strucObs2 = strucObs;
%IFAC_1DZE_extremeopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_1DZE_uinf11_noest_eopti_1em4/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4\workspace.mat')
for i = 1:Wp.sim.NN
    time21(i) = sol_array(i).score.CPUtime;
    RMSE21(i) = sol_array(i).score.RMSE_cline;
    maxError21(i) = sol_array(i).score.maxError;
    RMSE_flow21(i) = sol_array(i).score.RMSE_flow;
end
Wp21 = Wp; sol_array21 = sol_array; sys21 = sys; 
scriptOptions21 = scriptOptions; strucObs21 = strucObs;
%IFAC_2DZE_extremeopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf11_noest_eopti_1em4/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4\workspace.mat')
for i = 1:Wp.sim.NN
    time22(i) = sol_array(i).score.CPUtime;
    RMSE22(i) = sol_array(i).score.RMSE_cline;
    maxError22(i) = sol_array(i).score.maxError;
    RMSE_flow22(i) = sol_array(i).score.RMSE_flow;
end
Wp22 = Wp; sol_array22 = sol_array; sys22 = sys; 
scriptOptions22 = scriptOptions; strucObs22 = strucObs;

% x = [1:3];
% avg_time = [sum(time0)/length(time0),sum(time1)/length(time1),...
%             sum(time22)/length(time22)];
% RMSE_flow = [sum(RMSE_flow0)/length(RMSE_flow0),sum(RMSE_flow1)/length(RMSE_flow1),...
%             sum(RMSE_flow22)/length(RMSE_flow22)];
% maxError = [sum(maxError0)/length(maxError0),sum(maxError1/length(maxError1)),...
%             sum(maxError22)/length(maxError22)];
% figure, 
% yyaxis left, plot(x,avg_time,'--*'),
% ylabel('Average time (sec)')
% % yyaxis right, plot(x,RMSE_flow,':o')
% % ylabel('Average RMSE_flow')
% yyaxis right, plot(x,maxError,':o')
% ylabel('Average maxError')
% xticks(x)
% xticklabels({'OPEN','EnKF','2D','2Deopti1em4','3D','3Deopti1em4'})
% % xtickangle(90)
% title('Effect of estimation size') 
figure, plot(RMSE0), hold on, plot(RMSE1),plot(RMSE2), 
plot(RMSE21), plot(RMSE22)
legend('OPEN','EnKF:1D','EnKF:2D','DExKF_IFAC:1DE_extopti_1em4','DExKF_IFAC:2DE_extopti_1em4')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.1);
title('RMSE_cline'), xlabel('Time (sec)'), ylabel('RMSE')
figure, plot(maxError0), hold on, plot(maxError1),plot(maxError2),
plot(maxError21),plot(maxError22)
legend('OPEN','EnKF:1D','EnKF:2D','DExKF_IFAC:1DE_extopti_1em4','DExKF_IFAC:2DE_extopti_1em4')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.1);
title('maxError'), xlabel('Time (sec)'), ylabel('maxError')
figure, plot(RMSE_flow0), hold on, plot(RMSE_flow1),plot(RMSE_flow2),
plot(RMSE_flow21),plot(RMSE_flow22)
legend('OPEN','EnKF:1D','EnKF:2D','DExKF_IFAC:1DE_extopti_1em4','DExKF_IFAC:2DE_extopti_1em4')
lq = findobj(gcf,'type','line');
set(lq,'linewidth',1.1);
title('RMSE_flow'), xlabel('Time (sec)'), ylabel('RMSE_flow')

hold off
plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
plotWFObs( Wp2,sol_array2,sys2,scriptOptions2,strucObs2 );
plotWFObs( Wp21,sol_array21,sys21,scriptOptions21,strucObs21 );
plotWFObs( Wp22,sol_array22,sys22,scriptOptions22,strucObs22 );
%% Linearizing at different frequencies with 2D and U_inf = 5 and it is not estimated
clear all
close all
clc
%IFAC_2DZE_extopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    cRMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
    time1(i) = sol_array(i).score.CPUtime;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys; 
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%IFAC_2DZE_extopti_1em4_NL_10%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_ext_1em4_NL_10/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_ext_1em4_NL_10\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    cRMSE2(i) = sol_array(i).score.RMSE_cline;
    maxError2(i) = sol_array(i).score.maxError;
    RMSE_flow2(i) = sol_array(i).score.RMSE_flow;
    time2(i) = sol_array(i).score.CPUtime;
end
Wp2 = Wp; sol_array2 = sol_array; sys2 = sys; 
scriptOptions2 = scriptOptions; strucObs2 = strucObs;
%IFAC_2DZE_extopti_1em4_NL_20%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_ext_1em4_NL_20/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_ext_1em4_NL_20\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    cRMSE3(i) = sol_array(i).score.RMSE_cline;
    maxError3(i) = sol_array(i).score.maxError;
    RMSE_flow3(i) = sol_array(i).score.RMSE_flow;
    time3(i) = sol_array(i).score.CPUtime;
end
Wp3 = Wp; sol_array3 = sol_array; sys3 = sys; 
scriptOptions3 = scriptOptions; strucObs3 = strucObs;
%IFAC_2DZE_extopti_1em4_NL_50%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_ext_1em4_NL_50/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_ext_1em4_NL_50\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    cRMSE4(i) = sol_array(i).score.RMSE_cline;
    maxError4(i) = sol_array(i).score.maxError;
    RMSE_flow4(i) = sol_array(i).score.RMSE_flow;
    time4(i) = sol_array(i).score.CPUtime;
end
Wp4 = Wp; sol_array4 = sol_array; sys4 = sys; 
scriptOptions4 = scriptOptions; strucObs4 = strucObs;
%IFAC_2DZE_extopti_1em4_NL_100%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_ext_1em4_NL_100/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_ext_1em4_NL_100\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    cRMSE5(i) = sol_array(i).score.RMSE_cline;
    maxError5(i) = sol_array(i).score.maxError;
    RMSE_flow5(i) = sol_array(i).score.RMSE_flow;
    time5(i) = sol_array(i).score.CPUtime;
end
Wp5 = Wp; sol_array5 = sol_array; sys5 = sys; 
scriptOptions5 = scriptOptions; strucObs5 = strucObs;
figure, plot(cRMSE1), hold on, 
plot(cRMSE2), plot(cRMSE3), plot(cRMSE4), plot(cRMSE5), 
legend('DExKF_IFAC:2D_extopti_NL_inf','DExKF_IFAC:2D_extopti_NL_10','DExKF_IFAC:2D_extopti_NL_20','DExKF_IFAC:2D_extopti_NL_50','DExKF_IFAC:2D_extopti_NL_100')
title('RMSE_cline'), xlabel('Time (sec)'), ylabel('RMSE')
figure, plot(maxError1), hold on, 
plot(maxError2), plot(maxError3), plot(maxError4), plot(maxError5), 
legend('DExKF_IFAC:2D_extopti_NL_inf','DExKF_IFAC:2D_extopti_NL_10','DExKF_IFAC:2D_extopti_NL_20','DExKF_IFAC:2D_extopti_NL_50','DExKF_IFAC:2D_extopti_NL_100')
title('maxError'), xlabel('Time (sec)'), ylabel('maxError')
figure, plot(RMSE_flow1), hold on, 
plot(RMSE_flow2), plot(RMSE_flow3), plot(RMSE_flow4), plot(RMSE_flow5), 
legend('DExKF_IFAC:2D_extopti_NL_inf','DExKF_IFAC:2D_extopti_NL_10','DExKF_IFAC:2D_extopti_NL_20','DExKF_IFAC:2D_extopti_NL_50','DExKF_IFAC:2D_extopti_NL_100')
title('RMSE_flow'), xlabel('Time (sec)'), ylabel('RMSE_flow')
x = [1:5];
avg_time = [sum(time2),sum(time3),sum(time4),sum(time5),sum(time1)]./Wp.sim.NN;
RMSE_flow = [sum(RMSE_flow2),sum(RMSE_flow3),sum(RMSE_flow4),sum(RMSE_flow5),sum(RMSE_flow1)]./Wp.sim.NN;
maxError = [sum(maxError2),sum(maxError3),sum(maxError4),sum(maxError5),sum(maxError1)]./Wp.sim.NN;
figure, 
yyaxis left, plot(x,avg_time,'--*'),
ylabel('Average time (sec)')
% yyaxis right, plot(x,RMSE_flow,':o')
% ylabel('Average RMSE_flow')
yyaxis right, plot(x,maxError,':o')
ylabel('Average maxError')
xticks(x)
xticklabels({'Every 10 secs','Every 20 secs','Every 50 secs','Every 100 secs','1^{st} iteration'})
title('Effect of linearizing at different frequencies') 

% hold off
% plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
% plotWFObs( Wp2,sol_array2,sys2,scriptOptions2,strucObs2 );
% plotWFObs( Wp3,sol_array3,sys3,scriptOptions3,strucObs3 );
% plotWFObs( Wp4,sol_array4,sys4,scriptOptions4,strucObs4 );
% plotWFObs( Wp5,sol_array5,sys5,scriptOptions5,strucObs5 );
%% Different initial value for mixing length and it is not estimated (EnKF, DExKF_IFAC:1D,2D)
% With E
clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_sim_lmu_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    u0(i) = sol_array(i).site.u_Inf;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys; 
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_lmu_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u1(i) = sol_array(i).site.u_Inf;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys; 
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%IFAC_1DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_1DZE_lmu_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZE_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u2(i) = sol_array(i).site.u_Inf;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
    maxError2(i) = sol_array(i).score.maxError;
    RMSE_flow2(i) = sol_array(i).score.RMSE_flow;
end
Wp2 = Wp; sol_array2 = sol_array; sys2 = sys; 
scriptOptions2 = scriptOptions; strucObs2 = strucObs;
%IFAC_2DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_lmu_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim11,timeWFSim11,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u3(i) = sol_array(i).site.u_Inf;
    RMSE3(i) = sol_array(i).score.RMSE_cline;
    maxError3(i) = sol_array(i).score.maxError;
    RMSE_flow3(i) = sol_array(i).score.RMSE_flow;
end
Wp3 = Wp; sol_array3 = sol_array; sys3 = sys; 
scriptOptions3 = scriptOptions; strucObs3 = strucObs;
%IFAC_3DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_3DZE_lmu_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim12,timeWFSim12,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u4(i) = sol_array(i).site.u_Inf;
    RMSE4(i) = sol_array(i).score.RMSE_cline;
    maxError4(i) = sol_array(i).score.maxError;
    RMSE_flow4(i) = sol_array(i).score.RMSE_flow;
end
Wp4 = Wp; sol_array4 = sol_array; sys4 = sys; 
scriptOptions4 = scriptOptions; strucObs4 = strucObs;

figure, plot(RMSE0), hold on, 
plot(RMSE1),plot(RMSE2), plot(RMSE3), plot(RMSE4) 
legend('OPEN','EnKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:3D')
title('RMSE_cline'), xlabel('Time (sec)'), ylabel('RMSE')
figure, plot(maxError0), hold on, 
plot(maxError1),plot(maxError2), plot(maxError3), plot(maxError4)
legend('OPEN','EnKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:3D')
title('RMSE'), xlabel('Time (sec)'), ylabel('maxError')
figure, plot(RMSE_flow0), hold on, 
plot(RMSE_flow1),plot(RMSE_flow2), plot(RMSE_flow3), plot(RMSE_flow4)
legend('OPEN','EnKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:3D')
title('RMSE'), xlabel('Time (sec)'), ylabel('RMSE_flow')

hold off
plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
plotWFObs( Wp2,sol_array2,sys2,scriptOptions2,strucObs2 );
plotWFObs( Wp3,sol_array3,sys3,scriptOptions3,strucObs3 );
plotWFObs( Wp4,sol_array4,sys4,scriptOptions4,strucObs4 );
%% Different initial value for freestream velocity and it is not estimated (EnKF, DExKF_IFAC:1D,2D)
% Without E
clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_sim_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    u0(i) = sol_array(i).site.u_Inf;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys; 
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u1(i) = sol_array(i).site.u_Inf;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys; 
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_ExKF_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u2(i) = sol_array(i).site.u_Inf;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
    maxError2(i) = sol_array(i).score.maxError;
    RMSE_flow2(i) = sol_array(i).score.RMSE_flow;
end
Wp2 = Wp; sol_array2 = sol_array; sys2 = sys; 
scriptOptions2 = scriptOptions; strucObs2 = strucObs;
%IFAC_1DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_1DZ_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u3(i) = sol_array(i).site.u_Inf;
    RMSE3(i) = sol_array(i).score.RMSE_cline;
    maxError3(i) = sol_array(i).score.maxError;
    RMSE_flow3(i) = sol_array(i).score.RMSE_flow;
end
Wp3 = Wp; sol_array3 = sol_array; sys3 = sys; 
scriptOptions3 = scriptOptions; strucObs3 = strucObs;
%IFAC_2DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZ_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u4(i) = sol_array(i).site.u_Inf;
    RMSE4(i) = sol_array(i).score.RMSE_cline;
    maxError4(i) = sol_array(i).score.maxError;
    RMSE_flow4(i) = sol_array(i).score.RMSE_flow;
end
Wp4 = Wp; sol_array4 = sol_array; sys4 = sys; 
scriptOptions4 = scriptOptions; strucObs4 = strucObs;
%IFAC_4DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_4DZ_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_4DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u5(i) = sol_array(i).site.u_Inf;
    RMSE5(i) = sol_array(i).score.RMSE_cline;
    maxError5(i) = sol_array(i).score.maxError;
    RMSE_flow5(i) = sol_array(i).score.RMSE_flow;
end
Wp5 = Wp; sol_array5 = sol_array; sys5 = sys; 
scriptOptions5 = scriptOptions; strucObs5 = strucObs;
%CIN_4DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_CIN_4DZ_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_4DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim6,timeWFSim6,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u6(i) = sol_array(i).site.u_Inf;
    RMSE6(i) = sol_array(i).score.RMSE_cline;
    maxError6(i) = sol_array(i).score.maxError;
    RMSE_flow6(i) = sol_array(i).score.RMSE_flow;
end
Wp6 = Wp; sol_array6 = sol_array; sys6 = sys; 
scriptOptions6 = scriptOptions; strucObs6 = strucObs;
%CIN_1DC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_CIN_1DC_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_1DC_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim7,timeWFSim7,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u7(i) = sol_array(i).site.u_Inf;
    RMSE7(i) = sol_array(i).score.RMSE_cline;
    maxError7(i) = sol_array(i).score.maxError;
    RMSE_flow7(i) = sol_array(i).score.RMSE_flow;
end
Wp7 = Wp; sol_array7 = sol_array; sys7 = sys; 
scriptOptions7 = scriptOptions; strucObs7 = strucObs;
%CIN_2DC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_CIN_2DC_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_2DC_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim8,timeWFSim8,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u8(i) = sol_array(i).site.u_Inf;
    RMSE8(i) = sol_array(i).score.RMSE_cline;
    maxError8(i) = sol_array(i).score.maxError;
    RMSE_flow8(i) = sol_array(i).score.RMSE_flow;
end
Wp8 = Wp; sol_array8 = sol_array; sys8 = sys; 
scriptOptions8 = scriptOptions; strucObs8 = strucObs;
%CIN_1DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_CIN_1DZ_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_1DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim9,timeWFSim9,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u9(i) = sol_array(i).site.u_Inf;
    RMSE9(i) = sol_array(i).score.RMSE_cline;
    maxError9(i) = sol_array(i).score.maxError;
    RMSE_flow9(i) = sol_array(i).score.RMSE_flow;
end
Wp9 = Wp; sol_array9 = sol_array; sys9 = sys; 
scriptOptions9 = scriptOptions; strucObs9 = strucObs;
figure, plot(RMSE0), hold on, 
plot(RMSE1),plot(RMSE2), plot(RMSE3), plot(RMSE4), plot(RMSE5), 
plot(RMSE6), plot(RMSE7), plot(RMSE8), plot(RMSE9), 
legend('OPEN','EnKF','ExKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:4D','DExKF_CIN:4DZ','DExKF_CIN:1DC','DExKF_CIN:2DC','DExKF_CIN:1DZ')
title('RMSE_cline'), xlabel('Time (sec)'), ylabel('RMSE')
figure, plot(maxError0), hold on, 
plot(maxError1),plot(maxError2), plot(maxError3), plot(maxError4), plot(maxError5), 
plot(maxError6), plot(maxError7), plot(maxError8), plot(maxError9), 
legend('OPEN','EnKF','ExKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:4D','DExKF_CIN:4DZ','DExKF_CIN:1DC','DExKF_CIN:2DC','DExKF_CIN:1DZ')
title('RMSE'), xlabel('Time (sec)'), ylabel('maxError')
figure, plot(RMSE_flow0), hold on, 
plot(RMSE_flow1),plot(RMSE_flow2), plot(RMSE_flow3), plot(RMSE_flow4), plot(RMSE_flow5), 
plot(RMSE_flow6), plot(RMSE_flow7), plot(RMSE_flow8), plot(RMSE_flow9), 
legend('OPEN','EnKF','ExKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:4D','DExKF_CIN:4DZ','DExKF_CIN:1DC','DExKF_CIN:2DC','DExKF_CIN:1DZ')
title('RMSE'), xlabel('Time (sec)'), ylabel('RMSE_flow')

hold off
plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
plotWFObs( Wp2,sol_array2,sys2,scriptOptions2,strucObs2 );
plotWFObs( Wp3,sol_array3,sys3,scriptOptions3,strucObs3 );
plotWFObs( Wp4,sol_array4,sys4,scriptOptions4,strucObs4 );
plotWFObs( Wp5,sol_array5,sys5,scriptOptions5,strucObs5 );
plotWFObs( Wp6,sol_array6,sys6,scriptOptions6,strucObs6 );
plotWFObs( Wp7,sol_array7,sys7,scriptOptions7,strucObs7 );
plotWFObs( Wp8,sol_array8,sys8,scriptOptions8,strucObs8 );
plotWFObs( Wp9,sol_array9,sys9,scriptOptions9,strucObs9 );
%% Different initial value for freestream velocity and it is not estimated (EnKF, DExKF_IFAC:1D,2D)
% With E
clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_sim_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    time0(i) = sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys; 
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_EnKF_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time1(i) = sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys; 
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
% %ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf_noest/workspace.mat')
% % load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_ExKF_uinf_noest\workspace.mat')
% for i = 1:Wp.sim.NN
%     time4(i) = sol_array(i).score.CPUtime;
%     RMSE4(i) = sol_array(i).score.RMSE_cline;
%     maxError4(i) = sol_array(i).score.maxError;
%     RMSE_flow4(i) = sol_array(i).score.RMSE_flow;
% end
% Wp4 = Wp; sol_array4 = sol_array; sys4 = sys; 
% scriptOptions4 = scriptOptions; strucObs4 = strucObs;
% %IFAC_1DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_1DZE_uinf_noest/workspace.mat')
% % load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZE_uinf_noest\workspace.mat')
% % scriptOptions.powerForecast = 0;
% % [pwWFSim10,timeWFSim10,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     time10(i) = sol_array(i).score.CPUtime;
%     RMSE10(i) = sol_array(i).score.RMSE_cline;
%     maxError10(i) = sol_array(i).score.maxError;
%     RMSE_flow10(i) = sol_array(i).score.RMSE_flow;
% end
% Wp10 = Wp; sol_array10 = sol_array; sys10 = sys; 
% scriptOptions10 = scriptOptions; strucObs10 = strucObs;
%IFAC_2DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim11,timeWFSim11,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time11(i) = sol_array(i).score.CPUtime;
    RMSE11(i) = sol_array(i).score.RMSE_cline;
    maxError11(i) = sol_array(i).score.maxError;
    RMSE_flow11(i) = sol_array(i).score.RMSE_flow;
end
Wp11 = Wp; sol_array11 = sol_array; sys11 = sys; 
scriptOptions11 = scriptOptions; strucObs11 = strucObs;
%IFAC_3DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim12,timeWFSim12,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time12(i) = sol_array(i).score.CPUtime;
    RMSE12(i) = sol_array(i).score.RMSE_cline;
    maxError12(i) = sol_array(i).score.maxError;
    RMSE_flow12(i) = sol_array(i).score.RMSE_flow;
end
Wp12 = Wp; sol_array12 = sol_array; sys12 = sys; 
scriptOptions12 = scriptOptions; strucObs12 = strucObs;
% %NoFusion_1DCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_nofusion_1DC_uinf_noest\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim13,timeWFSim13,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u13(i) = sol_array(i).site.u_Inf;
%     RMSE13(i) = sol_array(i).score.RMSE_cline;
%     maxError13(i) = sol_array(i).score.maxError;
%     RMSE_flow13(i) = sol_array(i).score.RMSE_flow;
% end
% Wp13 = Wp; sol_array13 = sol_array; sys13 = sys; 
% scriptOptions13 = scriptOptions; strucObs13 = strucObs;
% %NoFusion_2DCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_nofusion_2DC_uinf_noest\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim14,timeWFSim14,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u14(i) = sol_array(i).site.u_Inf;
%     RMSE14(i) = sol_array(i).score.RMSE_cline;
%     maxError14(i) = sol_array(i).score.maxError;
%     RMSE_flow14(i) = sol_array(i).score.RMSE_flow;
% end
% Wp14 = Wp; sol_array14 = sol_array; sys14 = sys; 
% scriptOptions14 = scriptOptions; strucObs14 = strucObs;
%IFAC_1DZE_opti%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZE_uinf_noest_opti\workspace.mat')
% % scriptOptions.powerForecast = 0;
% % [pwWFSim15,timeWFSim15,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     time15(i) = sol_array(i).score.CPUtime;
%     RMSE15(i) = sol_array(i).score.RMSE_cline;
%     maxError15(i) = sol_array(i).score.maxError;
%     RMSE_flow15(i) = sol_array(i).score.RMSE_flow;
% end
% Wp15 = Wp; sol_array15 = sol_array; sys15 = sys; 
% scriptOptions15 = scriptOptions; strucObs15 = strucObs;
% %IFAC_2DZE_opti%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_opti\workspace.mat')
% % scriptOptions.powerForecast = 0;
% % [pwWFSim16,timeWFSim16,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     time16(i) = sol_array(i).score.CPUtime;
%     RMSE16(i) = sol_array(i).score.RMSE_cline;
%     maxError16(i) = sol_array(i).score.maxError;
%     RMSE_flow16(i) = sol_array(i).score.RMSE_flow;
% end
% Wp16 = Wp; sol_array16 = sol_array; sys16 = sys; 
% scriptOptions16 = scriptOptions; strucObs16 = strucObs;
% %IFAC_3DZE_opti%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_opti\workspace.mat')
% % scriptOptions.powerForecast = 0;
% % [pwWFSim17,timeWFSim17,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     time17(i) = sol_array(i).score.CPUtime;
%     RMSE17(i) = sol_array(i).score.RMSE_cline;
%     maxError17(i) = sol_array(i).score.maxError;
%     RMSE_flow17(i) = sol_array(i).score.RMSE_flow;
% end
% Wp17 = Wp; sol_array17 = sol_array; sys17 = sys; 
% scriptOptions17 = scriptOptions; strucObs17 = strucObs;
%IFAC_2DZE_superopti_1em2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em2/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em2\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim18,timeWFSim18,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time18(i) = sol_array(i).score.CPUtime;
    RMSE18(i) = sol_array(i).score.RMSE_cline;
    maxError18(i) = sol_array(i).score.maxError;
    RMSE_flow18(i) = sol_array(i).score.RMSE_flow;
end
Wp18 = Wp; sol_array18 = sol_array; sys18 = sys; 
scriptOptions18 = scriptOptions; strucObs18 = strucObs;
%IFAC_2DZE_superopti_1em3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em3/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em3\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim19,timeWFSim19,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time19(i) = sol_array(i).score.CPUtime;
    RMSE19(i) = sol_array(i).score.RMSE_cline;
    maxError19(i) = sol_array(i).score.maxError;
    RMSE_flow19(i) = sol_array(i).score.RMSE_flow;
end
Wp19 = Wp; sol_array19 = sol_array; sys19 = sys; 
scriptOptions19 = scriptOptions; strucObs19 = strucObs;
%IFAC_2DZE_superopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em4/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em4/workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim20,timeWFSim20,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time20(i) = sol_array(i).score.CPUtime;
    RMSE20(i) = sol_array(i).score.RMSE_cline;
    maxError20(i) = sol_array(i).score.maxError;
    RMSE_flow20(i) = sol_array(i).score.RMSE_flow;
end
Wp20 = Wp; sol_array20 = sol_array; sys20 = sys; 
scriptOptions20 = scriptOptions; strucObs20 = strucObs;
%IFAC_3DZE_superopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_superopti_1em4/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_superopti_1em4\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim21,timeWFSim21,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time21(i) = sol_array(i).score.CPUtime;
    RMSE21(i) = sol_array(i).score.RMSE_cline;
    maxError21(i) = sol_array(i).score.maxError;
    RMSE_flow21(i) = sol_array(i).score.RMSE_flow;
end
Wp21 = Wp; sol_array21 = sol_array; sys21 = sys; 
scriptOptions21 = scriptOptions; strucObs21 = strucObs;
x = [1:13];
avg_time = [sum(time0)/length(time0),sum(time1)/length(time1),...
    sum(time4)/length(time4),sum(time10)/length(time10),...
    sum(time15)/length(time11),sum(time11)/length(time15),...
    sum(time16)/length(time15),sum(time18)/length(time15),...
    sum(time19)/length(time15),sum(time20)/length(time15),...
    sum(time12)/length(time16),sum(time17)/length(time17),...
    sum(time21)/length(time15)];
RMSE_flow = [sum(RMSE_flow0)/length(RMSE_flow0),sum(RMSE_flow1)/length(RMSE_flow1),...
    sum(RMSE_flow4)/length(RMSE_flow4),sum(RMSE_flow10)/length(RMSE_flow10),...
    sum(RMSE_flow15)/length(RMSE_flow11),sum(RMSE_flow11)/length(RMSE_flow12),...
    sum(RMSE_flow16)/length(RMSE_flow15),sum(RMSE_flow18)/length(RMSE_flow15),...
    sum(RMSE_flow19)/length(RMSE_flow15),sum(RMSE_flow20)/length(RMSE_flow15),...
    sum(RMSE_flow12)/length(RMSE_flow16),sum(RMSE_flow17)/length(RMSE_flow17),...
    sum(RMSE_flow21)/length(RMSE_flow15)];
maxError = [sum(maxError0)/length(maxError0),sum(maxError1/length(maxError1)),...
    sum(maxError4)/length(maxError4),sum(maxError10)/length(maxError10),...
    sum(maxError15)/length(maxError11),sum(maxError11)/length(maxError12),...
    sum(maxError16)/length(maxError15),sum(maxError18)/length(maxError15),...
    sum(maxError19)/length(maxError15),sum(maxError20)/length(maxError15),...
    sum(maxError12)/length(maxError16),sum(maxError17)/length(maxError17),...
    sum(maxError21)/length(maxError15)];
figure, 
yyaxis left, plot(x,avg_time,'--*'),
ylabel('Average time (sec)')
% yyaxis right, plot(x,RMSE_flow,':o')
% ylabel('Average RMSE_flow')
yyaxis right, plot(x,maxError,':o')
ylabel('Average maxError')
xticks(x)
xticklabels({'OPEN','EnKF','ExKF','1D','1Dopti','2D','2Dopti','2Dsopti1em2','2Dsopti1em3','2Dsopti1em4','3D','3Dopti','3Dsopti1em4'})
xtickangle(90)
title('Effect of estimation size') 
figure,
plot(RMSE_flow0), hold on,plot(RMSE_flow1),plot(RMSE_flow4),
plot(RMSE_flow10),plot(RMSE_flow15),plot(RMSE_flow11),
plot(RMSE_flow16),plot(RMSE_flow18),plot(RMSE_flow19),
plot(RMSE_flow20),plot(RMSE_flow12),plot(RMSE_flow17),
plot(RMSE_flow21),
legend('OPEN','EnKF','ExKF','1D','1Dopti','2D','2Dopti','2Dsopti1em2','2Dsopti1em3','2Dsopti1em4','3D','3Dopti','3Dsopti1em4');
title('RMSE Flow'), xlabel('time (sec)'),ylabel('RMSE')
% hold off
% plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
% plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
% plotWFObs( Wp4,sol_array4,sys4,scriptOptions4,strucObs4 );
% plotWFObs( Wp10,sol_array10,sys10,scriptOptions10,strucObs10 );
% plotWFObs( Wp11,sol_array11,sys11,scriptOptions11,strucObs11 );
% plotWFObs( Wp12,sol_array12,sys12,scriptOptions12,strucObs12 );
% plotWFObs( Wp15,sol_array15,sys15,scriptOptions15,strucObs15 );
% plotWFObs( Wp16,sol_array16,sys16,scriptOptions16,strucObs16 );
% plotWFObs( Wp17,sol_array17,sys17,scriptOptions17,strucObs17 );
% plotWFObs( Wp18,sol_array18,sys18,scriptOptions18,strucObs18 );
% plotWFObs( Wp19,sol_array19,sys19,scriptOptions19,strucObs19 );
% plotWFObs( Wp20,sol_array20,sys20,scriptOptions20,strucObs20 );
% plotWFObs( Wp21,sol_array21,sys21,scriptOptions21,strucObs21 );
%IFAC_2DZE_extremeopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_extremeopti_1em4\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim22,timeWFSim22,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time22(i) = sol_array(i).score.CPUtime;
    RMSE22(i) = sol_array(i).score.RMSE_cline;
    maxError22(i) = sol_array(i).score.maxError;
    RMSE_flow22(i) = sol_array(i).score.RMSE_flow;
end
Wp22 = Wp; sol_array22 = sol_array; sys22 = sys; 
scriptOptions22 = scriptOptions; strucObs22 = strucObs;
%IFAC_3DZE_extremeopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_extremeopti_1em4/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_extremeopti_1em4\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim23,timeWFSim23,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    time23(i) = sol_array(i).score.CPUtime;
    RMSE23(i) = sol_array(i).score.RMSE_cline;
    maxError23(i) = sol_array(i).score.maxError;
    RMSE_flow23(i) = sol_array(i).score.RMSE_flow;
end
Wp23 = Wp; sol_array23 = sol_array; sys23 = sys; 
scriptOptions23 = scriptOptions; strucObs23 = strucObs;
%IFAC_2DZE_extremeopti_1em4_Q1em1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_eopti_1em4_Q1em1/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_eopti_1em4_Q1em1\workspace.mat')
scriptOptions.powerForecast = 0;
for i = 1:Wp.sim.NN
    time24(i) = sol_array(i).score.CPUtime;
    RMSE24(i) = sol_array(i).score.RMSE_cline;
    maxError24(i) = sol_array(i).score.maxError;
    RMSE_flow24(i) = sol_array(i).score.RMSE_flow;
end
Wp24 = Wp; sol_array24 = sol_array; sys24 = sys; 
scriptOptions24 = scriptOptions; strucObs24 = strucObs;
%IFAC_3DZE_extremeopti_1em4_Q1em1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_eopti_1em4_Q1em1/workspace.mat')
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_eopti_1em4_Q1em1\workspace.mat')
scriptOptions.powerForecast = 0;
for i = 1:Wp.sim.NN
    time25(i) = sol_array(i).score.CPUtime;
    RMSE25(i) = sol_array(i).score.RMSE_cline;
    maxError25(i) = sol_array(i).score.maxError;
    RMSE_flow25(i) = sol_array(i).score.RMSE_flow;
end
Wp25 = Wp; sol_array25 = sol_array; sys25 = sys; 
scriptOptions25 = scriptOptions; strucObs25 = strucObs;

x = [1:6];
avg_time = [sum(time0)/length(time0),sum(time1)/length(time1),...
    sum(time11)/length(time15),sum(time22)/length(time15),...
    sum(time12)/length(time16),sum(time23)/length(time17)];
RMSE_flow = [sum(RMSE_flow0)/length(RMSE_flow0),sum(RMSE_flow1)/length(RMSE_flow1),...
    sum(RMSE_flow11)/length(RMSE_flow4),sum(RMSE_flow12)/length(RMSE_flow10),...
    sum(RMSE_flow12)/length(RMSE_flow11),sum(RMSE_flow23)/length(RMSE_flow12)];
maxError = [sum(maxError0)/length(maxError0),sum(maxError1/length(maxError1)),...
    sum(maxError11)/length(maxError11),sum(maxError22)/length(maxError22),...
    sum(maxError12)/length(maxError12),sum(maxError23)/length(maxError23)];
figure, 
yyaxis left, plot(x,avg_time,'--*'),
ylabel('Average time (sec)')
% yyaxis right, plot(x,RMSE_flow,':o')
% ylabel('Average RMSE_flow')
yyaxis right, plot(x,maxError,':o')
ylabel('Average maxError')
xticks(x)
xticklabels({'OPEN','EnKF','2D','2Deopti1em4','3D','3Deopti1em4'})
% xtickangle(90)
title('Effect of estimation size') 
figure,
plot(RMSE_flow0), hold on,plot(RMSE_flow1),plot(RMSE_flow11),
plot(RMSE_flow22),plot(RMSE_flow12),plot(RMSE_flow23),
legend('OPEN','EnKF','2D','2Deopti1em4','3D','3Deopti1em4');
title('RMSE Flow'), xlabel('time (sec)'),ylabel('RMSE')
hold off
plotWFObs( Wp22,sol_array22,sys22,scriptOptions22,strucObs22 );
plotWFObs( Wp23,sol_array23,sys23,scriptOptions23,strucObs23 );

figure, plot(RMSE0), hold on, plot(RMSE1), %plot(RMSE4), plot(RMSE10),
plot(RMSE11),plot(RMSE12),
% plot(RMSE15),plot(RMSE16),plot(RMSE17)
plot(RMSE18),plot(RMSE19),plot(RMSE20),plot(RMSE21)
plot(RMSE22),plot(RMSE23)
legend('OPEN','EnKF','DExKF_IFAC:2DE','DExKF_IFAC:3DE','DExKF_IFAC:2DE_superopti_1em2','DExKF_IFAC:2DE_superopti_1em3','DExKF_IFAC:2DE_superopti_1em4','DExKF_IFAC:3DE_superopti_1em4','DExKF_IFAC:2DE_extopti_1em4','DExKF_IFAC:3DE_extopti_1em4')
title('RMSE_cline'), xlabel('Time (sec)'), ylabel('RMSE')
figure, plot(maxError0), hold on, plot(maxError1), %plot(maxError4), plot(maxError10),
plot(maxError11),plot(maxError12), 
% plot(maxError15),plot(maxError16),plot(maxError17)
plot(maxError18),plot(maxError19),plot(maxError20),plot(maxError21)
plot(maxError22),plot(maxError23)
legend('OPEN','EnKF','DExKF_IFAC:2DE','DExKF_IFAC:3DE','DExKF_IFAC:2DE_superopti_1em2','DExKF_IFAC:2DE_superopti_1em3','DExKF_IFAC:2DE_superopti_1em4','DExKF_IFAC:3DE_superopti_1em4','DExKF_IFAC:2DE_extopti_1em4','DExKF_IFAC:3DE_extopti_1em4')
title('maxError'), xlabel('Time (sec)'), ylabel('maxError')
figure, plot(RMSE_flow0), hold on, plot(RMSE_flow1), %plot(RMSE_flow4),plot(RMSE_flow10),
plot(RMSE_flow11),plot(RMSE_flow12), 
% plot(RMSE_flow15),plot(RMSE_flow16),plot(RMSE_flow17)
plot(RMSE_flow18),plot(RMSE_flow19),plot(RMSE_flow20),plot(RMSE_flow21)
plot(RMSE_flow22),plot(RMSE_flow23)
legend('OPEN','EnKF','DExKF_IFAC:2DE','DExKF_IFAC:3DE','DExKF_IFAC:2DE_superopti_1em2','DExKF_IFAC:2DE_superopti_1em3','DExKF_IFAC:2DE_superopti_1em4','DExKF_IFAC:3DE_superopti_1em4','DExKF_IFAC:2DE_extopti_1em4','DExKF_IFAC:3DE_extopti_1em4')
title('RMSE_flow'), xlabel('Time (sec)'), ylabel('RMSE_flow')

hold off
plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
% % plotWFObs( Wp4,sol_array4,sys4,scriptOptions4,strucObs4 );
% % plotWFObs( Wp10,sol_array10,sys10,scriptOptions10,strucObs10 );
plotWFObs( Wp11,sol_array11,sys11,scriptOptions11,strucObs11 );
plotWFObs( Wp12,sol_array12,sys12,scriptOptions12,strucObs12 );
% % plotWFObs( Wp15,sol_array15,sys15,scriptOptions15,strucObs15 );
% % plotWFObs( Wp16,sol_array16,sys16,scriptOptions16,strucObs16 );
% % plotWFObs( Wp17,sol_array17,sys17,scriptOptions17,strucObs17 );
plotWFObs( Wp18,sol_array18,sys18,scriptOptions18,strucObs18 );
plotWFObs( Wp19,sol_array19,sys19,scriptOptions19,strucObs19 );
plotWFObs( Wp20,sol_array20,sys20,scriptOptions20,strucObs20 );
plotWFObs( Wp21,sol_array21,sys21,scriptOptions21,strucObs21 );
plotWFObs( Wp22,sol_array22,sys22,scriptOptions22,strucObs22 );
plotWFObs( Wp23,sol_array23,sys23,scriptOptions23,strucObs23 );
plotWFObs( Wp24,sol_array24,sys24,scriptOptions24,strucObs24 );
plotWFObs( Wp25,sol_array25,sys25,scriptOptions25,strucObs25 );
%% Comparing the effect of different estimation size (EnKF, DExKF_IFAC:1D,2D,3D,4D)
clear all
close all
clc

disp('--------------------------')
disp('EnKF')
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline0 = RMSE/(Wp.sim.NN - 1)
avg_maxError0 = maxError/(Wp.sim.NN - 1)
avg_time0 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline1 = RMSE/(Wp.sim.NN - 1)
avg_maxError1 = maxError/(Wp.sim.NN - 1)
avg_time1 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:2D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline2 = RMSE/(Wp.sim.NN - 1)
avg_maxError2 = maxError/(Wp.sim.NN - 1)
avg_time2 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:3D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE3(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline3 = RMSE/(Wp.sim.NN - 1)
avg_maxError3 = maxError/(Wp.sim.NN - 1)
avg_time3 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:4D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_4DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE4(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline4 = RMSE/(Wp.sim.NN - 1)
avg_maxError4 = maxError/(Wp.sim.NN - 1)
avg_time4 = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
plot(RMSE3),
plot(RMSE4),
legend('EnKF','DExKF:1D','DExKF:2D','DExKF:3D','DExKF:4D')
D = Wp.turbine.Drotor;
dia1    = D;
dia2    = 2*D;
dia3    = 3*D;
dia4    = 4*D;
dia     = [dia1,dia2,dia3,dia4];
avg_time= [avg_time1,avg_time2,avg_time3,avg_time4];
avg_RMSE_cline = [avg_RMSE_cline1,avg_RMSE_cline2,avg_RMSE_cline3,avg_RMSE_cline4];
figure, 
yyaxis left, plot(dia,avg_time,'-.*')
ylabel('Average time (sec)')
yyaxis right,plot(dia,avg_RMSE_cline,':o')
axis([D-(D/10) 4*D+(D/10), 0.3 0.5])
ylabel('Average RMSE Cline (sec)');
xlabel('Radius (m)')
xticks(dia)
xticklabels({'1D','2D','3D','4D'})
title('Effect of different estimation size')

%% Comparing different fusing algorithms (EnKF, DExKF_1D:IFAC,CIN,CI,EI,ICI)
clear all
close all
clc

disp('--------------------------')
disp('EnKF')
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline0 = RMSE/(Wp.sim.NN - 1)
avg_maxError0 = maxError/(Wp.sim.NN - 1)
avg_time0 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline1 = RMSE/(Wp.sim.NN - 1)
avg_maxError1 = maxError/(Wp.sim.NN - 1)
avg_time1 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_CIN')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_1DZ\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline2 = RMSE/(Wp.sim.NN - 1)
avg_maxError2 = maxError/(Wp.sim.NN - 1)
avg_time2 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% disp('--------------------------')
% disp('DExKF:1D_CI')
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CI_1DZ\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% RMSE        = 0;
% maxError    = 0;
% time        = 0;
% for i = 1:Wp.sim.NN
%     RMSE        = RMSE + sol_array(i).score.RMSE_cline;
%     maxError    = maxError + sol_array(i).score.maxError;
%     time        = time + sol_array(i).score.CPUtime;
%     RMSE3(i) = sol_array(i).score.RMSE_cline;
% end
% avg_RMSE_cline3 = RMSE/(Wp.sim.NN - 1)
% avg_maxError3 = maxError/(Wp.sim.NN - 1)
% avg_time3 = time/(Wp.sim.NN - 1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% disp('--------------------------')
% disp('DExKF:1D_EI')
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_EI_1DZ\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% RMSE        = 0;
% maxError    = 0;
% time        = 0;
% for i = 1:Wp.sim.NN
%     RMSE        = RMSE + sol_array(i).score.RMSE_cline;
%     maxError    = maxError + sol_array(i).score.maxError;
%     time        = time + sol_array(i).score.CPUtime;
%     RMSE4(i) = sol_array(i).score.RMSE_cline;
% end
% avg_RMSE_cline4 = RMSE/(Wp.sim.NN - 1)
% avg_maxError4 = maxError/(Wp.sim.NN - 1)
% avg_time4 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_ICI')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_ICI_1DZ\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE5(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline5 = RMSE/(Wp.sim.NN - 1)
avg_maxError5 = maxError/(Wp.sim.NN - 1)
avg_time5 = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
% plot(RMSE3),
% plot(RMSE4),
plot(RMSE5),
legend('EnKF','DExKF:IFAC','DExKF:CIN','DExKF:ICI')
% legend('EnKF','DExKF:IFAC','DExKF:CIN','DExKF:CI','DExKF:EI','DExKF:ICI')
x = [1,2,3,4];
avg_time= [avg_time0,avg_time1,avg_time2,avg_time5];
avg_RMSE_cline = [avg_RMSE_cline0,avg_RMSE_cline1,avg_RMSE_cline2,avg_RMSE_cline5];
figure, 
yyaxis left, plot(x,avg_time,'-.*')
ylabel('Average time (sec)')
yyaxis right,plot(x,avg_RMSE_cline,':o')
% axis([D-(D/10) 4*D+(D/10), 0.3 0.5])
ylabel('Average RMSE Cline (sec)');
xlabel('Type of Filter')
xticks(x)
xticklabels({'EnKF','DExKF: IFAC','DExKF: CI(Naive)','DExKF: ICI'})
title('Effect of different fusion algorithms')

%% Comparing the effects of linearizing the non-linear system at different frequencies (DExKF:1D_IFAC: 10,20,50,100,Inf)
clear all
close all
clc

disp('--------------------------')
disp('EnKF')
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline0 = RMSE/(Wp.sim.NN - 1)
avg_maxError0 = maxError/(Wp.sim.NN - 1)
avg_time0 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: 10')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_10\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline1 = RMSE/(Wp.sim.NN - 1)
avg_maxError1 = maxError/(Wp.sim.NN - 1)
avg_time1 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: 20')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_20\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline2 = RMSE/(Wp.sim.NN - 1)
avg_maxError2 = maxError/(Wp.sim.NN - 1)
avg_time2 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: 50')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_50\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE3(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline3 = RMSE/(Wp.sim.NN - 1)
avg_maxError3 = maxError/(Wp.sim.NN - 1)
avg_time3 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: 100')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_100\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE4(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline4 = RMSE/(Wp.sim.NN - 1)
avg_maxError4 = maxError/(Wp.sim.NN - 1)
avg_time4 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: Inf')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE5(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline5 = RMSE/(Wp.sim.NN - 1)
avg_maxError5 = maxError/(Wp.sim.NN - 1)
avg_time5 = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
plot(RMSE3),
plot(RMSE4),
plot(RMSE5),
legend('EnKF','DExKF:10','DExKF:20','DExKF:50','DExKF:100','DExKF:Inf')
x = [1,2,3,4,5,6];
avg_time= [avg_time0,avg_time1,avg_time2,avg_time3,avg_time4,avg_time5];
avg_RMSE_cline = [avg_RMSE_cline0,avg_RMSE_cline1,avg_RMSE_cline2,avg_RMSE_cline3,avg_RMSE_cline4,avg_RMSE_cline5];
figure, 
yyaxis left, plot(x,avg_time,'-.*')
ylabel('Average time (sec)')
yyaxis right,plot(x,avg_RMSE_cline,':o')
% axis([D-(D/10) 4*D+(D/10), 0.3 0.5])
ylabel('Average RMSE Cline (sec)');
xlabel('Type of Filter')
xticks(x)
xticklabels({'EnKF','DExKF:10','DExKF:20','DExKF:50','DExKF:100','DExKF:Inf'})
title('Effect of linearizing the system at different time-steps')

%% Disturbances added to the initial flow field (EnKF, DExKF_1D:IFAC,ICI)
clear all
close all
clc

load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE3(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZ_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE4(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZ_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE5(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim6,timeWFSim6,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE6(i) = sol_array(i).score.RMSE_cline;
end
figure, plot(RMSE1), hold on, plot(RMSE2), plot(RMSE3), plot(RMSE4), plot(RMSE5), plot(RMSE6)
legend('EnKF','DExKF_IFAC:1D','DExKF: ICI','DExKF_IFAC:2D','DExKF_IFAC:3D','sim')

%% Disturbances added to the output (EnKF, DExKF:1D_IFAC)
clear all
close all
clc

load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_noise_obs_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
% plot(u), hold on, legend('EnKF', 'DExKF')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_noise_obs_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
figure, plot(RMSE1), hold on, plot(RMSE2)
legend('EnKF','DExKF')
%% Estimate the freestream flow (EnKF, DExKF:1D_IFAC)
clear all
close all
clc

load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_uinf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
u1 = [];
for i = 1:Wp.sim.NN
    u1(i) = sol_array(i).site.u_Inf;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
% plot(u), hold on, legend('EnKF', 'DExKF')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_uinf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
u2 = [];
for i = 1:Wp.sim.NN
    u2(i) = sol_array(i).site.u_Inf;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
figure, plot(u1), hold on, 
plot(u2), plot(Wp.site.actual_u_Inf*ones(1,Wp.sim.NN)) 
legend('EnKF','DExKF'), axis([0 2000 0 9])
title('Estimation of Freestream Velcoity'), xlabel('Time (sec)'), ylabel('U_Inf')
figure, plot(RMSE1), hold on, plot(RMSE2)
legend('EnKF','DExKF')
title('RMSE'), xlabel('Time (sec)'), ylabel('RMSE')
%%
clear all
close all
clc
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_ExKF/ExKF_2turb_workspace.mat')
% load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_ExKF\ExKF_2turb_workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim1,timeWFSim1,pwLES,timeLES] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_EnKF/EnKF_2turb_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_EnKF\EnKF_2turb_workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 2/dexkf_2turb_IFAC_1D_workspace.mat')
% load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 2\dexkf_2turb_IFAC_1D_workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 3/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 3\dexkf_2turb_IFAC_1D_workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 4/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 4\dexkf_2turb_IFAC_1D_workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC2DZ/trail 2/dexkf_2turb_IFAC_2D_workspace.mat')
% load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC2DZ\trail 2\dexkf_2turb_IFAC_2D_workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim6,timeWFSim6,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_dexkf_ICI1DZ\trail 2\dexkf_2turb_ICI_1D_workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim7,timeWFSim7,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
hold off

% e11 = pwWFSim1(1,:) - pwLES(1,:);
% e12 = pwWFSim1(2,:) - pwLES(2,:);
% 
% e21 = pwWFSim2(1,:) - pwLES(1,:);
% e22 = pwWFSim2(2,:) - pwLES(2,:);
% 
% e31 = pwWFSim3(1,:) - pwLES(1,:);
% e32 = pwWFSim3(2,:) - pwLES(2,:);
% 
% e41 = pwWFSim4(1,:) - pwLES(1,:);
% e42 = pwWFSim4(2,:) - pwLES(2,:);
% 
% e51 = pwWFSim5(1,:) - pwLES(1,:);
% e52 = pwWFSim5(2,:) - pwLES(2,:);
% 
% figure,plot(timeWFSim1,e11), hold on, plot(timeWFSim1,e21),plot(timeWFSim1,e31), plot(timeWFSim1,e41), hold off
% figure,plot(timeWFSim1,e12), hold on, plot(timeWFSim1,e22),plot(timeWFSim1,e32), plot(timeWFSim1,e42), hold off
% 
% ee11 = sum(e11)/length(e11)
% ee21 = sum(e21)/length(e21)
% ee31 = sum(e31)/length(e31)
% ee41 = sum(e41)/length(e41)
% 
% ee12 = sum(e12)/length(e12)
% ee22 = sum(e22)/length(e22)
% ee32 = sum(e32)/length(e32)
% ee42 = sum(e42)/length(e42)

%% 2Turbine

clear all
close all
clc

disp('--------------------------')
disp('ExKF: Linearize every iterations')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_ExKF/ExKF_2turb_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_ExKF\ExKF_2turb_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('EnKF')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_EnKF/EnKF_2turb_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_EnKF\EnKF_2turb_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 1D + IFAC: Trail 1 (Unoptimized + linearize once every 5 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 1/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 1\dexkf_2turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 1D + IFAC: Trail 2 (Optimized + linearize once every 50 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 2/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 2\dexkf_2turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 1D + IFAC: Trail 3 (Optimized + linearize once every 5 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 3/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 3\dexkf_2turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 1D + IFAC: Trail 4 (Optimized + linearize at k = 1 + Conv. KF)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 3/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 4\dexkf_2turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 2D + IFAC')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC2DZ/trail 1/dexkf_2turb_IFAC_2D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC2DZ\trail 1\dexkf_2turb_IFAC_2D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 2D + IFAC: Trail 2 (Optimized + linearize once every 50 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC2DZ/trail 2/dexkf_2turb_IFAC_2D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC2DZ\trail 2\dexkf_2turb_IFAC_2D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 3D + IFAC')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC3DZ/dexkf_2turb_IFAC_3D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC3DZ\dexkf_2turb_IFAC_3D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 4D + IFAC')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC4DZ/dexkf_2turb_IFAC_4D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC4DZ\dexkf_2turb_IFAC_4D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 1D + ICI: Trail 1 (Unoptimized + linearize once every 5 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_dexkf_ICI1DZ/dexkf_2turb_ICI_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_dexkf_ICI1DZ\trail 1\dexkf_2turb_ICI_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 1D + ICI: Trail 2 (Optimized + linearize at k = 1)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_dexkf_ICI1DZ/dexkf_2turb_ICI_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_dexkf_ICI1DZ\trail 2\dexkf_2turb_ICI_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%% 9Turbine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('9TURBINE: SIM')
clear all
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/9Turbine/apc_9turb_alm_turb/apc_9turb_alm_turb_sim/9turb_sim_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('9TURBINE: DExKF: 1D + IFAC')
clear all
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/9Turbine/apc_9turb_alm_turb/apc_9turb_alm_turb_dexkf_IFAC1DZ/Trail 1/dexkf_9turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

disp('--------------------------')
disp('9TURBINE: DExKF: 1D + IFAC')
clear all
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/9Turbine/apc_9turb_alm_turb/apc_9turb_alm_turb_dexkf_IFAC1DZ/Trail 2/dexkf_9turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

disp('--------------------------')
disp('9TURBINE (apc_9turb_adm_noturb): DExKF: 1D + IFAC')
clear all
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/9Turbine/apc_9turb_adm_noturb/apc_9turb_adm_noturb_dexkf_IFAC1DZ/workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );