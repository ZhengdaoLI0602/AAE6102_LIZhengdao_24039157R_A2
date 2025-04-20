
%% Task 2 
clearvars
load('Outputs/urban_matlab.mat')
%%%==============Load Skymask=================%%%
skymask = readmatrix("FGI-GSRx-2.0.1\AAE6102_ZhengdaoLI_A2\skymask_A1_urban.csv");
angles     = skymask(2:end,1);
elevations = skymask(2:end,2);

angles_smooth = linspace(min(angles), max(angles), 1e4); % 1e4 points for smoothness
smoothing_param = 0.99; % 0 (no smoothing) to 1 (max smoothing)
pp = csaps(angles, elevations, smoothing_param);
elevations_smooth = fnval(pp, angles_smooth);
settings.sys.skymask = [angles_smooth', elevations_smooth'];
%%%==============Load Skymask=================%%%

% Execute navigation
settings.sys.ls_type = 0; % 0 for OLS;  1 for WLS
settings.sys.raim_type = 2; % 0 for raim;  1 for weighted raim; 2 for no raim
settings.sys.skymask_type = 0; % 0 for no skymask;  1 for skymask filtering
[obsData, satData, navData] = doNavigation(obsData, settings, ephData);






%% Task 3 
clearvars
load('Outputs/opensky_matlab.mat')
% Execute navigation
settings.sys.ls_type = 1; % 0 for OLS;  1 for WLS
settings.sys.raim_type = 1; % 0 for raim;  1 for weighted raim; 2 for no raim
[obsData, satData, navData] = doNavigation(obsData, settings, ephData);







%% Written by Zhengdao LI (zhengda0.li@connect.polyu.hk)
close all

% test_ekf
%% Save and load data
% save urban_nav_OLS navData settings
% save urban_nav_OLS_skymask navData settings
% save nav_WLS navData settings
% save nav_EKF navData 

clearvars
% load("nav_OLS.mat")
load("urban_nav_OLS.mat")
solu_ols = navData;
load("urban_nav_OLS_skymask.mat")
solu_ols_skymask = navData;
% load("nav_WLS.mat")
% solu_wls = navData;

NrEpoch    = size(navData, 2);
GT_llh     = [settings.nav.trueLat, settings.nav.trueLong, settings.nav.trueHeight];
GT_ecef    = llh2ecef(GT_llh.*[pi/180, pi/180, 1]);
R_mat      = R_ecef_enu(GT_llh.*[pi/180, pi/180, 1]); % transform from LLH to ENU coordinate

pos_llh = (1:3);
vel_enu = (4:6);
pos_enu = (7:9);
pos_raim = (10:13);

for i = 1: NrEpoch % loop through each epoch
    % OLS: (pos_llh, vel_enu)
    SOLU.ols(i, pos_llh) = solu_ols{1,i}.Pos.LLA;
    SOLU.ols(i, vel_enu) = (R_mat * solu_ols{1,i}.Vel.xyz')'; % ((3,3) * (1,3)')' = (1,3)
    SOLU.ols(i, pos_enu) = (R_mat * (solu_ols{1,i}.Pos.xyz - GT_ecef)' )';
    % SOLU.ols(i, pos_raim) = solu_ols{1,i}.Pos.raim;

    % OLS with skymask: (pos_llh, vel_enu)
    SOLU.ols_skymask(i, pos_llh) = solu_ols_skymask{1,i}.Pos.LLA;
    SOLU.ols_skymask(i, vel_enu) = (R_mat * solu_ols_skymask{1,i}.Vel.xyz')'; % ((3,3) * (1,3)')' = (1,3)
    SOLU.ols_skymask(i, pos_enu) = (R_mat * (solu_ols_skymask{1,i}.Pos.xyz - GT_ecef)' )';
    % SOLU.ols_skymask(i, pos_raim) = solu_ols_skymask{1,i}.Pos.raim;

    % % WLS: (pos_llh, vel_enu)
    % SOLU.wls(i, pos_llh) = solu_wls{1,i}.Pos.LLA;
    % SOLU.wls(i, vel_enu) = (R_mat * solu_wls{1,i}.Vel.xyz')';
    % SOLU.wls(i, pos_enu) = (R_mat * (solu_wls{1,i}.Pos.xyz - GT_ecef)' )';
    % SOLU.wls(i, pos_raim) = solu_wls{1,i}.Pos.raim;
end

% Location
% figure; hold on; grid on
% scatter(GT_llh(2), GT_llh(1), 150, 'filled', 'g', "pentagram",'MarkerEdgeColor', 'k', 'DisplayName', 'Grount Truth')
% scatter(SOLU.ols(:, pos_llh(2)), SOLU.ols(:, pos_llh(1)), 30, 'filled', 'r', 'DisplayName', 'OLS')
% % scatter(SOLU.wls(:, pos_llh(2)), SOLU.wls(:, pos_llh(1)), 30, 'filled', 'b', 'DisplayName', 'WLS')
% scatter(SOLU.ols_skymask(:, pos_llh(2)), SOLU.ols_skymask(:, pos_llh(1)), 30, 'filled', 'b', 'DisplayName', 'OLS with Skymask')
% legend('show');

figure; 
geobasemap('satellite');
geoscatter(SOLU.ols(:, pos_llh(1)), SOLU.ols(:, pos_llh(2)), 30, 'filled', 'r', 'DisplayName', 'OLS'); hold on
geoscatter(SOLU.ols_skymask(:, pos_llh(1)), SOLU.ols_skymask(:, pos_llh(2)), 30, 'filled', 'b', 'DisplayName', 'OLS with Skymask'); hold on
geoscatter(GT_llh(1), GT_llh(2), 300, 'filled', 'g', "pentagram",'MarkerEdgeColor', 'k', 'DisplayName', 'Grount Truth'); hold on
legend('show');



err_3d_ols = vecnorm(SOLU.ols(:, pos_enu),2,2);
err_3d_ols_skymask = vecnorm(SOLU.ols_skymask(:, pos_enu),2,2);

mean(err_3d_ols)
std(err_3d_ols)
mean(err_3d_ols_skymask)
std(err_3d_ols_skymask)

SOLU.ols_skymask(:, pos_enu)


fde_info = SOLU.ols(:, pos_raim);
nav_info = SOLU.ols;

figure; 
subplot(1,3,1)
hold on
plot(fde_info(:, 1), 'b', LineWidth=1); % sqrt WSSE
plot(fde_info(:, 2), 'r', LineWidth=1); % Threshold
xlabel("Time (Epoch)")
ylabel("Test Statistics (m)")
title("Detection results for WLS")

subplot(1,3,2)
plot(fde_info(:, 3), 'g', LineWidth=1); % PL
xlabel("Time (Epoch)")
ylabel("PL (m)")
title("PL for WLS")

err_3d = vecnorm(nav_info(:, pos_enu), 2, 2);
pl_3d = fde_info(:, 3);
subplot(1,3,3)
hold on;
% Scatter plot of Protection Level vs. Position Error
scatter(err_3d, pl_3d, 'filled', 'b');
ylabel('Protection Level (m)');
xlabel('Position Error (m)');
% title('Stanford Plot for OLS');
grid on;
% Add alert limit lines
yline(50, '--black', 'LineWidth', 1.5);
xline(50, '--black', 'LineWidth', 1.5);
% Add 1:1 reference line (ideal scenario)
plot([0, 60], [0, 60], 'k--', 'LineWidth', 1);
title("Stanford chart for WLS")


% Position Error
figure; 
subplot(3,1,1)
hold on; grid on
plot(1:NrEpoch, SOLU.ols(:,pos_enu(1)), 'b-*', 'DisplayName', 'E', LineWidth=1);
plot(1:NrEpoch, SOLU.ols(:,pos_enu(2)), 'g-+', 'DisplayName', 'N', LineWidth=1); 
plot(1:NrEpoch, SOLU.ols(:,pos_enu(3)), 'r-o', 'DisplayName', 'U', LineWidth=1);  
legend('show');
ylabel("Error (m)")
title("OLS")

subplot(3,1,2)
hold on; grid on
plot(1:NrEpoch, SOLU.wls(:,pos_enu(1)), 'b-*', 'DisplayName', 'E', LineWidth=1);
plot(1:NrEpoch, SOLU.wls(:,pos_enu(2)), 'g-+', 'DisplayName', 'N', LineWidth=1); 
plot(1:NrEpoch, SOLU.wls(:,pos_enu(3)), 'r-o', 'DisplayName', 'U', LineWidth=1);  
legend('show');
ylabel("Error (m)")
title("WLS")

subplot(3,1,3)
hold on; grid on
plot(1:NrEpoch, SOLU.ekf(:,pos_enu(1)), 'b-*', 'DisplayName', 'E', LineWidth=1);
plot(1:NrEpoch, SOLU.ekf(:,pos_enu(2)), 'g-+', 'DisplayName', 'N', LineWidth=1); 
plot(1:NrEpoch, SOLU.ekf(:,pos_enu(3)), 'r-o', 'DisplayName', 'U', LineWidth=1);  
legend('show');
ylabel("Error (m)")
title("EKF")
xlabel("Time (Epoch)")





% Velocity
figure; 
subplot(3,1,1)
hold on; grid on
plot(1:NrEpoch, SOLU.ols(:,vel_enu(1)), 'b-*', 'DisplayName', 'E', LineWidth=1);
plot(1:NrEpoch, SOLU.ols(:,vel_enu(2)), 'g-+', 'DisplayName', 'N', LineWidth=1); 
plot(1:NrEpoch, SOLU.ols(:,vel_enu(3)), 'r-o', 'DisplayName', 'U', LineWidth=1);  
legend('show');
ylabel("Velocity (m/s)")
title("OLS")

subplot(3,1,2)
hold on; grid on
plot(1:NrEpoch, SOLU.wls(:,vel_enu(1)), 'b-*', 'DisplayName', 'E', LineWidth=1);
plot(1:NrEpoch, SOLU.wls(:,vel_enu(2)), 'g-+', 'DisplayName', 'N', LineWidth=1); 
plot(1:NrEpoch, SOLU.wls(:,vel_enu(3)), 'r-o', 'DisplayName', 'U', LineWidth=1);  
legend('show');
ylabel("Velocity (m/s)")
title("WLS")

subplot(3,1,3)
hold on; grid on
plot(1:NrEpoch, SOLU.ekf(:,vel_enu(1)), 'b-*', 'DisplayName', 'E', LineWidth=1);
plot(1:NrEpoch, SOLU.ekf(:,vel_enu(2)), 'g-+', 'DisplayName', 'N', LineWidth=1); 
plot(1:NrEpoch, SOLU.ekf(:,vel_enu(3)), 'r-o', 'DisplayName', 'U', LineWidth=1);  
legend('show');
ylabel("Velocity (m/s)")
title("EKF")

xlabel("Time (Epoch)")

