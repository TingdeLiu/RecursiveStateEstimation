%% Recursive State Estimation Exercise SoSe 2024 %%
% group number:  4                                %
% participant 1: Zhuoyue Xu                       %
% participant 2: Shaohao Yang                     %
% participant 3: Tingde Liu                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;clc;close all;format long g

%% loading the data
addpath(genpath(pwd),'data')
addpath(pwd,'filters')
addpath(pwd,'functions')
addpath(pwd,'results');
savePath = [pwd,'/results'];

load('GroundTruth_X_Y_ConstraintID.mat');   % ground truth
load('Observations_d1_a1_d2_a2.mat');   % measurements

%% store the true states in a matrix called "trueStates" with [x y ID] -> [m m -] %% 
trueStates = GroundTruth_X_Y_ConstraintID;

%% store the measurements in a matrix called "l" with [d1 a1 d2 a2] -> [m rad m rad] %%
l = Observations_d1_a1_d2_a2;

clear GroundTruth_X_Y_ConstraintID
clear Observations_d1_a1_d2_a2

%% given information in the data-sheet %%
TS1 = [135.54;98.79];   % position of TS1
TS2 = [110.00;90.00];   % positiion of TS2
R = 9.00;   % radious of a sector of the trajectory
sigma.stdx = 1e-4; % [m]
sigma.stdy = sigma.stdx;
sigma.vx = 1e-1;   % [m/s]
sigma.vy = sigma.vx;
sigma.w = 1e-2; % noise value 
delta_t = 1;    % [s]
ep_first = 1;
ep_last = size(l,1);

%% form the initial state vector and call it "Xinit" and cofactor matrix and call it "Qxxinit" %%   
Xinit = [trueStates(1,1);trueStates(1,2);0;0];   % [x,y,vx,vy]
v = [sigma.stdx,sigma.stdy,sigma.vx,sigma.vy];
Qxxinit = diag(v.^2);

%% form the VCM of the process noise and call it Qww %%
v = [sigma.w,sigma.w,sigma.w,sigma.w];
Qww = diag(v.^2);

%% form the transition matrix and call it "Phi" %%
Phi = [1,0,delta_t,0;...
0,1,0,delta_t;...
0,0,1,0;...
0,0,0,1];

%% derive the trajectory from the observations %%
trajT1 = [TS1(1,1) + l(:,1).*cos(l(:,2)) TS1(2,1) + l(:,1).*sin(l(:,2))];    % [X-coordinates  Y-coordinates] by means of TS1
    % trajT2 -> [X-coordinates  Y-coordinates] by means of TS2 
trajT2 = [TS2(1,1) + l(:,3).*cos(l(:,4)) TS2(2,1) + l(:,3).*sin(l(:,4))];
%% average of the total statian observations as the measured trajectory of the vehicle %%
T_Mean_X = mean([trajT1(:,1) trajT2(:,1)],2);
T_Mean_Y = mean([trajT1(:,2) trajT2(:,2)],2);
T_Mean = [T_Mean_X T_Mean_Y];

%% Filters %%
%% Extended Kalman Filter (EKF) %%
[EKF_finalStates,EKF_finalQxx] = EKF(TS1,TS2,l,Xinit,Qxxinit,Qww,Phi,trueStates,false,ep_first,ep_last);

%% Extended Kalman Filter with State Constraints (EKFC) %%
[EKFC_finalStates,EKFC_finalQxx] = EKF(TS1,TS2,l,Xinit,Qxxinit,Qww,Phi,trueStates,true,ep_first,ep_last);

%% Iterated Extended Kalman FIlter (IEKF) %%
noIter = 10; % number of iterations 
[IEKF_finalStates,IEKF_finalQxx] = IEKF(TS1,TS2,l,Xinit,Qxxinit,Qww,Phi,trueStates,false,noIter,ep_first,ep_last);

%% Iterated Extended Kalman FIlter with State Constraints (IEKFC) %%
noIter = 10; % number of iterations 
[IEKFC_finalStates,IEKFC_finalQxx] = IEKF(TS1,TS2,l,Xinit,Qxxinit,Qww,Phi,trueStates,true,noIter,ep_first,ep_last);

%% Unscented Kalman Filter (UKF) %%
gamma = 1;  % scale factor
[UKF_finalStates,UKF_finalQxx] = UKF(TS1,TS2,l,Xinit,Qxxinit,Qww,Phi,gamma,ep_first,ep_last);

%% Ensemble Kalman Filter (EnKF) %%
numberEnsemble = 1000;  % number of ensembles to be generated
[EnKF_finalStates,EnKF_finalQxx] = EnKF(TS1,TS2,l,Xinit,Qxxinit,Qww,Phi,numberEnsemble,ep_first,ep_last);

%% Particle Filter (PF) %%
noParticle = 10000;
[PF_finalStates,PF_finalQxx] = PF(TS1,TS2,l,Xinit,Qxxinit,Qww,Phi,noParticle,ep_first,ep_last);

%% RMSE %%
% form a matrix (matDiffFilt2True_EKF) that shows the differences between
% the estimated states from EKF and the true states
matDiffFilt2True_EKF = sqrt(sum((EKF_finalStates(:,1:2)-trueStates(:,1:2)).^2,2)); % sqrt(delta_x^2+delta_y^2)

% form a matrix (matDiffFilt2True_EKFC) that shows the differences between
% the estimated states from EKFC and the true states
matDiffFilt2True_EKFC = sqrt(sum((EKFC_finalStates(:,1:2)-trueStates(:,1:2)).^2,2));

% form a matrix (matDiffFilt2True_IEKF) that shows the differences between
% the estimated states from IEKF and the true states
matDiffFilt2True_IEKF = sqrt(sum((IEKF_finalStates(:,1:2)-trueStates(:,1:2)).^2,2));

% form a matrix (matDiffFilt2True_IEKFC) that shows the differences between
% the estimated states from IEKFC and the true states
matDiffFilt2True_IEKFC = sqrt(sum((IEKFC_finalStates(:,1:2)-trueStates(:,1:2)).^2,2));

% form a matrix (matDiffFilt2True_UKF) that shows the differences between
% the estimated states from UKF and the true states
matDiffFilt2True_UKF = sqrt(sum((UKF_finalStates(:,1:2)-trueStates(:,1:2)).^2,2));

% form a matrix (matDiffFilt2True_EnKF) that shows the differences between
% the estimated states from EnkF and the true states
matDiffFilt2True_EnKF = sqrt(sum((EnKF_finalStates(:,1:2)-trueStates(:,1:2)).^2,2));

% form a matrix (matDiffFilt2True_PF) that shows the differences between
% the estimated states from PF and the true states
matDiffFilt2True_PF = sqrt(sum((PF_finalStates(:,1:2)-trueStates(:,1:2)).^2,2));

% form a matrix (matDiffFilt2True_TS) that shows the differences between
% the average of the total statian observations (T_mean) and the true
% states
matDiffFilt2True_TS = sqrt(sum((T_Mean-trueStates(:,1:2)).^2,2));

% calculate the cumulative RMSE for the EKF solutions (call it "RMSE_EKF"),
% the EKFC solutions (call it "RMSE_EKFC"), the IEKF 
% solutions (call it "RMSE_IEKF"), the IEKFC solutions 
% (call it "RMSE_IEKFC"), the UKF solutions (call it 
% "RMSE_UKF"), the EnKF solutions (call it "RMSE_EnKF"), the PF solutions 
% (call it "RMSE_PF") and the total station measurements 
% (call it "RMSE_TS")
for epoch = ep_first:ep_last
    if epoch == ep_first
        RMSE_EKF(epoch,:) = sqrt(matDiffFilt2True_EKF(1:epoch,1).^2 / epoch);
        RMSE_TS(epoch,:) = sqrt(matDiffFilt2True_TS(1:epoch,1).^2 / epoch);
        RMSE_EKFC(epoch,:) = sqrt(matDiffFilt2True_EKFC(1:epoch,1).^2 / epoch);
        RMSE_IEKF(epoch,:) = sqrt(matDiffFilt2True_IEKF(1:epoch,1).^2 / epoch);
        RMSE_IEKFC(epoch,:) = sqrt(matDiffFilt2True_IEKFC(1:epoch,1).^2 / epoch);
        RMSE_UKF(epoch,:) = sqrt(matDiffFilt2True_UKF(1:epoch,1).^2 / epoch);
        RMSE_EnKF(epoch,:) = sqrt(matDiffFilt2True_EnKF(1:epoch,1).^2 / epoch);
        RMSE_PF(epoch,:) = sqrt(matDiffFilt2True_PF(1:epoch,1).^2 / epoch);
    else
        RMSE_EKF(epoch,:) = sqrt( sum( matDiffFilt2True_EKF(1:epoch,1).^2) / epoch);
        RMSE_TS(epoch,:) = sqrt( sum( matDiffFilt2True_TS(1:epoch,1).^2) / epoch);
        RMSE_EKFC(epoch,:) = sqrt( sum( matDiffFilt2True_EKFC(1:epoch,1).^2) / epoch);
        RMSE_IEKF(epoch,:) = sqrt(sum( matDiffFilt2True_IEKF(1:epoch,1).^2)/ epoch);
        RMSE_IEKFC(epoch,:) = sqrt(sum( matDiffFilt2True_IEKFC(1:epoch,1).^2)/ epoch);
        RMSE_UKF(epoch,:) = sqrt(sum( matDiffFilt2True_UKF(1:epoch,1).^2)/ epoch);
        RMSE_EnKF(epoch,:) = sqrt(sum( matDiffFilt2True_EnKF(1:epoch,1).^2)/ epoch);
        RMSE_PF(epoch,:) = sqrt(sum( matDiffFilt2True_PF(1:epoch,1).^2)/ epoch);
    end
end

% take the RMSE value in the last epoch and scale it to 1000 (multiply by
% 1000) -> for the EKF call it "endRMSE_EKF" for the EKF solutions, 
% "endRMSE_EKFC" for the EKFC solutions, "endRMSE_IEKF" for 
% the IEKF solutions, "endRMSE_IEKFC" for the IEKFC solutions, 
% "endRMSE_UKF" for the UKF solutions, "endRMSE_EnKF" for the EnKF 
% solutions, "endRMSE_PF" for the PF solutions, and "endRMSE_TS" for the 
% total station measurements
endRMSE_EKF = 1000 * sqrt( sum( matDiffFilt2True_EKF.^2) / numel(matDiffFilt2True_EKF(:,1)))'; 
endRMSE_TS = 1000 * sqrt( sum( matDiffFilt2True_TS.^2) / numel(matDiffFilt2True_TS(:,1)))';
endRMSE_EKFC = 1000 * sqrt( sum( matDiffFilt2True_EKFC.^2) / numel(matDiffFilt2True_EKFC(:,1)))'; 
endRMSE_IEKF = 1000 * sqrt( sum( matDiffFilt2True_IEKF.^2) / numel(matDiffFilt2True_IEKF(:,1)))'; 
endRMSE_IEKFC = 1000 * sqrt( sum( matDiffFilt2True_IEKFC.^2) / numel(matDiffFilt2True_IEKFC(:,1)))'; 
endRMSE_UKF = 1000 * sqrt( sum( matDiffFilt2True_UKF.^2) / numel(matDiffFilt2True_UKF(:,1)))'; 
endRMSE_EnKF = 1000 * sqrt( sum( matDiffFilt2True_EnKF.^2) / numel(matDiffFilt2True_EnKF(:,1)))'; 
endRMSE_PF = 1000 * sqrt( sum( matDiffFilt2True_PF.^2) / numel(matDiffFilt2True_PF(:,1)))'; 
%% Plots %%
epochs = (ep_first:ep_last)';

figure()
hold all
grid on
plot(T_Mean(:,1),T_Mean(:,2),'*g')
plot(trueStates(:,1),trueStates(:,2),'*r')
plot(TS1(1,1),TS1(2,1),'om','MarkerSize',10,'LineWidth',2)
plot(TS2(1,1),TS2(2,1),'oc','Markersize',10,'LineWidth',2)
plot(EKF_finalStates(:,1),EKF_finalStates(:,2),'*m')
plot(EKFC_finalStates(:,1),EKFC_finalStates(:,2),'*c')
plot(IEKF_finalStates(:,1),IEKF_finalStates(:,2),'*y')
plot(IEKFC_finalStates(:,1),IEKFC_finalStates(:,2),'*k')
plot(UKF_finalStates(:,1),UKF_finalStates(:,2),'*','color',[127/255,0/255,255/255])
plot(EnKF_finalStates(:,1),EnKF_finalStates(:,2),'*','color',[204/255,0/255,102/255])
plot(PF_finalStates(:,1),PF_finalStates(:,2),'*','color',[96/255,96/255,96/255])
xlim([100 145])
ylim([80 115])
title('taken trajectory of the vehicle')
xlabel('X [m]')
ylabel('Y [m]')
legend('from observations','ground truth','TS1','TS2','EKF','EKFC','IEKF','IEKFC','UKF','EnKF','PF')
set(gca,'FontSize',20,'FontWeight','bold')
hold off
figName = "trajectories.fig";
savefig(fullfile(savePath,figName));
close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold all
grid on
plot(epochs,RMSE_EKF*1000,'m','LineWidth',2)
plot(epochs,RMSE_EKFC*1000,'c','LineWidth',2)
plot(epochs,RMSE_IEKF*1000,'y','LineWidth',2)
plot(epochs,RMSE_IEKFC*1000,'k','LineWidth',2)
plot(epochs,RMSE_UKF*1000,'color',[127/255,0/255,255/255],'LineWidth',2)
plot(epochs,RMSE_EnKF*1000,'color',[204/255,0/255,102/255],'LineWidth',2)
plot(epochs,RMSE_PF*1000,'color',[96/255,96/255,96/255],'LineWidth',2)
plot(epochs,RMSE_TS*1000,'g','LineWidth',2)
xlabel('epochs [s]')
ylabel('RMSE [mm]')
legend('EKF','EKFC','IEKF','IEKFC','UKF','EnKF','PF','TS')
set(gca,'FontSize',20,'FontWeight','bold')
hold off
figName = "RMSE.fig";
savefig(fullfile(savePath,figName));
close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
names = categorical({'EKF','EKFC','IEKF','IEKFC','UKF','EnKF','PF','TS'});
names = reordercats(names,{'EKF','EKFC','IEKF','IEKFC','UKF','EnKF','PF','TS'});
RMSE = [endRMSE_EKF endRMSE_EKFC endRMSE_IEKF endRMSE_IEKFC endRMSE_UKF endRMSE_EnKF endRMSE_PF endRMSE_TS];
bar(names,RMSE)
text(1:length(RMSE),RMSE,num2str(RMSE',3),'vert','bottom','horiz','center'); 
title('cumulative RMSE in the last epoch')
ylabel('RMSE [mm]')
set(gca,'FontSize',20,'FontWeight','bold')
figName = "barRMSE";
savefig(fullfile(savePath,figName));
close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general configuration
genConfCheck.true = trueStates;
genConfCheck.obs = l;
genConfCheck.initX = Xinit;
genConfCheck.initQxx = Qxxinit;
genConfCheck.Qww = Qww;
genConfCheck.Phi = Phi;
genConfCheck.trajT1 = trajT1;
genConfCheck.trajT2 = trajT2;

% EKF
EKFCheck.X = EKF_finalStates;
EKFCheck.Qxx = EKF_finalQxx;

% EKFC 
EKFCCheck.X = EKFC_finalStates;
EKFCCheck.Qxx = EKFC_finalQxx;

% IEKF
IEKFCheck.X = IEKF_finalStates;
IEKFCheck.Qxx = IEKF_finalQxx;

% IEKFC
IEKFCCheck.X = IEKFC_finalStates;
IEKFCCheck.Qxx = IEKFC_finalQxx;

% UKF 
UKFCheck.X = UKF_finalStates;
UKFCheck.Qxx = UKF_finalQxx;

% EnKF
EnKFCheck.X = EnKF_finalStates;
EnKFCheck.Qxx = EnKF_finalQxx;

% PF
PFCCheck.X = PF_finalStates;
PFCCheck.Qxx = PF_finalQxx;

% RMSE
rootMeanCheck.diffEKF = matDiffFilt2True_EKF;
rootMeanCheck.diffEKFC = matDiffFilt2True_EKFC;
rootMeanCheck.diffIEKF = matDiffFilt2True_IEKF;
rootMeanCheck.diffIEKFC = matDiffFilt2True_IEKFC;
rootMeanCheck.diffUKF = matDiffFilt2True_UKF;
rootMeanCheck.diffEnKF = matDiffFilt2True_EnKF;
rootMeanCheck.diffPF = matDiffFilt2True_PF;
rootMeanCheck.diffTS = matDiffFilt2True_TS;
rootMeanCheck.RMSEEKF = RMSE_EKF;
rootMeanCheck.RMSEEKFC = RMSE_EKFC;
rootMeanCheck.RMSEIEKF = RMSE_IEKF;
rootMeanCheck.RMSEIEKFC = RMSE_IEKFC;
rootMeanCheck.RMSEUKF = RMSE_UKF;
rootMeanCheck.RMSEEnKF = RMSE_EnKF;
rootMeanCheck.RMSEPF = RMSE_PF;
rootMeanCheck.RMSETS = RMSE_TS;

% saving
checkSolutions.general = genConfCheck;
checkSolutions.filterEKF = EKFCheck;
checkSolutions.filterEKFC = EKFCCheck;
checkSolutions.filterIEKF = IEKFCheck;
checkSolutions.filterIEKFC = IEKFCCheck;
checkSolutions.filterUKF = UKFCheck;
checkSolutions.filterEnKF = EnKFCheck;
checkSolutions.filterPF = PFCCheck;
checkSolutions.eval = rootMeanCheck;
save(fullfile(savePath,'checkSolutions.mat'),'checkSolutions')