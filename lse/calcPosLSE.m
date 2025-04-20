%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright 2015-2021 Finnish Geospatial Research Institute FGI, National
%% Land Survey of Finland. This file is part of FGI-GSRx software-defined
%% receiver. FGI-GSRx is a free software: you can redistribute it and/or
%% modify it under the terms of the GNU General Public License as published
%% by the Free Software Foundation, either version 3 of the License, or any
%% later version. FGI-GSRx software receiver is distributed in the hope
%% that it will be useful, but WITHOUT ANY WARRANTY, without even the
%% implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%% See the GNU General Public License for more details. You should have
%% received a copy of the GNU General Public License along with FGI-GSRx
%% software-defined receiver. If not, please visit the following website 
%% for further information: https://www.gnu.org/licenses/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Pos, collect_fault] = calcPosLSE(obs, sat, allSettings, Pos, collect_fault)
function [Pos, obs] = calcPosLSE(obs, sat, allSettings, Pos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function calculates the Least Square Solution
%
% Inputs:
%   obs             - Observations for one epoch
%   sat             - Satellite positions and velocities for one epoch
%   allSettings     - receiver settings
%   pos             - Initial position for the LSE 
%
% Outputs:
%   Pos             - receiver position and receiver clock error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% test on Logistic noise assumption
% Define gradient and Hessian
t = @(X,y,h) y-h*X;   %(1*1) - (1*N) * (N*1)
ft = @(X,y,h) log(exp(t(X,y,h))+exp(-t(X,y,h)));
% grad_ft = @(X,y,h) (exp(t(X,y,h))-exp(-t(X,y,h))) ./ (exp(t(X,y,h))+exp(-t(X,y,h)));
% hess_ft = @(X,y,h) 4 ./ (exp(t(X,y,h))+exp(-t(X,y,h)))^2;
grad_ft = @(t) (exp(t)-exp(-t)) ./ (exp(t)+exp(-t));
hess_ft = @(t) 4 ./ (exp(t)+exp(-t))^2;

% % Newton's method
% X_est = zeros(n, 1); % Initial guess
% tol = 1e-6;
% max_iter = 50;
% 
% 
% for iter = 1:max_iter
%     grad = 0;
%     hess = 0;
%     for sv_idx = 1: size(sv_matrix,1)
%         grad = grad + grad_ft(X_est, y(sv_idx), sv_matrix(sv_idx, :));  
%         hess = hess + hess_ft(X_est, y(sv_idx), sv_matrix(sv_idx, :)); 
%     end
%     delta_X = -hess \ grad; % Newton step
% 
%     X_est = X_est + delta_X;
%     if norm(delta_X) < tol
%         fprintf('Converged at iteration %d\n', iter);
%         break;
%     end
% end
%%

% Set starting point
pos = Pos.xyz;
Pos.bValid = false;

% Constants
scal_factor = 1;
WGS84oe = allSettings.const.EARTH_WGS84_ROT;
SPEED_OF_LIGHT = allSettings.const.SPEED_OF_LIGHT/scal_factor;

% Temporary variables
rcvr_clock_corr = 0;
nmbOfIterations = 10;

% Total number of signals enabled
nrOfSignals = allSettings.sys.nrOfSignals;


idx_to_isolate = -1;
whether_compute_PL = false;  
whether_abort_operation = false;



while true
    % Init clock elements in pos vector
    pos(4:3+nrOfSignals) = zeros; 
    % Iteratively find receiver position 
    for iter = 1:nmbOfIterations
        ind = 0;     
        % Loop over all signals
        for signalNr = 1:allSettings.sys.nrOfSignals
            
            % Extract signal acronym
            signal = allSettings.sys.enabledSignals{signalNr};        
            
            % Loop over all channels
            for channelNr = 1:obs.(signal).nrObs
                if(obs.(signal).channel(channelNr).bObsOk) && (ind+1 ~= idx_to_isolate)
                    ind = ind + 1; % Index for valid obervations     

                    % Calculate Azimuth angle and ELevation angle
                    ecef_sv       = sat.(signal).channel(channelNr).Pos;
                    llh_sv        = ecef2llh(ecef_sv);
                    llh_usr_true  = [allSettings.nav.trueLat, allSettings.nav.trueLong, allSettings.nav.trueHeight]...
                                    .* [pi/180, pi/180, 1];
                    lat_u  = llh_usr_true(1);   
                    lon_u = llh_usr_true(2);
                    lat_sv = llh_sv(1);  
                    lon_sv = llh_sv(2);
                    tan_AZ(ind) = sin(lon_sv-lon_u) / (cos(lat_u)*tan(lat_sv) - sin(lat_u)*cos(lon_sv-lon_u));

                    ecef_usr_true = llh2ecef(llh_usr_true);
                    enu_sv = ecef2enu(ecef_sv, ecef_usr_true);
                    sin_EL(ind) = enu_sv(3)/norm(enu_sv);

                    if allSettings.sys.skymask_type == 0
                        obs.(signal).channel(channelNr).skyMaskOk = true;
                    else
                        [~, min_i] = min(abs( atan(tan_AZ(ind))*180/pi-allSettings.sys.skymask(:,1) ));
                        if asin(sin_EL(ind))*180/pi < allSettings.sys.skymask(min_i,2)
                            obs.(signal).channel(channelNr).skyMaskOk = false;
                            ind = ind - 1;
                            continue;
                        end
                    end
                          
                    % Pseudorrange of satellites
                    pseudo_range(ind) = obs.(signal).channel(channelNr).corrP / scal_factor;
                                
                    % Calculate range to satellite
                    % dx=sat.(signal).channel(channelNr).Pos(1) / scal_factor -pos(1); % commented by me
                    % dy=sat.(signal).channel(channelNr).Pos(2) / scal_factor -pos(2);
                    % dz=sat.(signal).channel(channelNr).Pos(3) / scal_factor -pos(3);   
                    ecef_sv = ecef_sv / scal_factor;
                    dx=ecef_sv(1) - pos(1);
                    dy=ecef_sv(2) - pos(2);
                    dz=ecef_sv(3) - pos(3);  
                    range(ind)=sqrt(dx^2+dy^2+dz^2); % This is the calculated range to the satellites
    
                    % Direction cosines
                    sv_matrix(ind,1) = dx/range(ind);
                    sv_matrix(ind,2) = dy/range(ind);
                    sv_matrix(ind,3) = dz/range(ind);
                    sv_matrix(ind,3+signalNr) = 1;
    
    
                    % Calculate Carrier-to-Noise Ratio
                    CN0(ind) = mean(obs.(signal).channel(channelNr).CN0(2:end));
                    
                    % Total clock correction term (m). */
                    %clock_correction = c*(sv_pos.dDeltaTime - eph(info.PRN).group_delay);
                    clock_correction = 0;
                    
                    % First compute the SV's earth rotation correction
                    rhox = ecef_sv(1) - pos(1);
                    rhoy = ecef_sv(2) - pos(2);
                    % EarthRotCorr(ind) = WGS84oe / SPEED_OF_LIGHT * (sat.(signal).channel(channelNr).Pos(2)/scal_factor*rhox - sat.(signal).channel(channelNr).Pos(1)/scal_factor*rhoy);
                    EarthRotCorr(ind) = WGS84oe / SPEED_OF_LIGHT * (ecef_sv(2)*rhox - ecef_sv(1)*rhoy);

                    % Total propagation delay.
                    propagation_delay(ind) = range(ind) + EarthRotCorr(ind);  %propagation_delay(ind) = range(ind) + EarthRotCorr(ind) - clock_correction;
    
                    % Correct the pseudoranges also (because we corrected rcvr stamp)
                    % pseudo_range(ind)  = pseudo_range(ind) - SPEED_OF_LIGHT * rcvr_clock_corr; % commented by myself
                    omp.dRange(ind) = pseudo_range(ind) - propagation_delay(ind);
                    Res(ind) = omp.dRange(ind) - pos(3 + signalNr) * SPEED_OF_LIGHT;                
                end
            end
            nrSatsUsed(signalNr) = ind;
        end

        if size(sv_matrix,1) < 4
            break
            whether_abort_operation = true;
        end
    
        % Set Weighting matrix for WLS
        clear W;
        a = 30; A = 50; F = 20; T = 50;
        if allSettings.sys.ls_type == 1 % WLS
            for ind = 1: nrSatsUsed
                if CN0(ind) < T
                    W_diag(ind) = 1/sin_EL(ind)^2*( 10^(-(CN0(ind)-T)/a)*((A/10^(-(F-T)/a)-1)*(CN0(ind)-T)/(F-T)+1) );
                else
                    W_diag(ind) = 1;
                end
            end
            W = diag(1 ./ W_diag);
        else           % OLS
            W = 1;
        end
      
        W = W ./ sum(sum(W));
    
        % This is the actual solutions to the LSE optimisation problem
        clear H;
        clear dR; 
        H=sv_matrix;%(1:5,:);
        % dR=omp.dRange;%(1:5);
        dR=Res;
        DeltaPos = (H'*W*H)^(-1)*H'*W*dR'; % DeltaPos=(H'*H)^(-1)*H'*dR';

        % if norm(DeltaPos(1:3)) < 0.5
        %     %% Logistic noise assumption (2025.04.07)
        %     % clear y;
        %     % grad = zeros(3+nrOfSignals,1);
        %     % hess = zeros(3+nrOfSignals, 3+nrOfSignals);
        %     % y = dR';
        %     % for sv_idx = 1: size(H,1)
        %     %     t = y(sv_idx) - H(sv_idx, :) * DeltaPos;
        %     %     % t = t/1e6;
        %     %     % grad = grad + grad_ft(pos', y(sv_idx), H(sv_idx, :));  
        %     %     % hess = hess + hess_ft(pos', y(sv_idx), H(sv_idx, :)); 
        %     %     if abs(t) <= 709
        %     %         grad_ft = (exp(t)-exp(-t)) ./ (exp(t)+exp(-t));
        %     %         hess_ft = 4 ./ (exp(t)+exp(-t))^2;
        %     %     else
        %     %         if t > 0
        %     %             grad_ft = 1;
        %     %         else
        %     %             grad_ft = -1;
        %     %         end
        %     %         hess_ft = 0;
        %     %     end
        %     %     grad = grad + grad_ft * (-H(sv_idx, :)');  
        %     %     hess = hess + hess_ft * (H(sv_idx, :)' * H(sv_idx, :)); 
        %     % end
        %     % DeltaPos = -hess \ grad; % Newton step (t is large, hess=0, then no update)
        % 
        %      %% Gauss-Newton assumption
        %     clear y;
        %     y = dR';
        % 
        %     % Initialize residual vector
        %     J = zeros(size(H,1), size(H,2));
        %     t = zeros(size(H,1), 1);
        % 
        %     for sv_idx = 1: size(H,1)
        %         t(sv_idx) = y(sv_idx) - H(sv_idx, :) * DeltaPos;
        % 
        %         if abs(t(sv_idx)) <= 709
        %             grad_ft = (exp(t(sv_idx))-exp(-t(sv_idx))) ./ (exp(t(sv_idx))+exp(-t(sv_idx)));
        %         else
        %             if t(sv_idx) > 0
        %                 grad_ft = 1;
        %             else
        %                 grad_ft = -1;
        %             end
        %         end
        % 
        %         % Compute Jacobian for this residual
        %         J(sv_idx, :) = grad_ft * (-H(sv_idx, :));
        %     end
        %     if abs(t) > 709*ones(size(H,1), 1)
        %         gi = t;
        %     else
        %         gi = log(exp(t)+exp(-t));
        %     end
        %     DeltaPos = (J'*J)^(-1)*J' * gi;  % the correction seems not correct for user clock bias
        %     % Update gradient and Hessian approximations
        % end
        % DeltaPos
        %
    
        % Updating the position with the solution
        pos(1)=pos(1)-DeltaPos(1);
        pos(2)=pos(2)-DeltaPos(2);
        pos(3)=pos(3)-DeltaPos(3);
        

        % Update the clock offsets for all systems
        pos(4:end) = DeltaPos(4:end) / SPEED_OF_LIGHT; % In seconds
    end    
    
    
    % Weighted RAIM (2025.04.05)
    if allSettings.sys.raim_type == 2 % no RAIM: end the while loop
        break
    end
    Nr_sat = size(H,1);
    [Updated_mat, Detect_results] = chi2_detector(dR', W, H, Nr_sat, allSettings);
    if Detect_results.fault_confirmed == 0
        whether_compute_PL = true; 
        % no fault detected
        break
    else
        all_fault_check = 1e-6 * ones(Nr_sat,1);
        for sv_idx = 1: Nr_sat
            isolation_mat = ones(Nr_sat,1);
            isolation_mat(sv_idx) = 0;  % check idx^th sv 
            [~, Detect_test] = chi2_detector(dR', W, H, Nr_sat, allSettings, isolation_mat);
            all_fault_check(sv_idx) = Detect_test.fault_confirmed ;
        end
        if sum(all_fault_check) == Nr_sat-1 % having one and only one fault detected [1,1,1,0,1,1]
            idx_to_isolate = find(all_fault_check==0);
            whether_compute_PL = true;
            % step to the WLS a second time without isolated satellite
            continue
        else
            % detected more than one faulty measurements, cannot used for positioning 
            whether_abort_operation = true;
            break
        end
    end
end


% if whether_compute_PL
%     PL = compute_PL(Updated_mat, Detect_results);
% else
%     PL = NaN;
% end
% 
% collect_fault = [collect_fault; [Detect_results.WSSE_sqrt, Detect_results.Thres, PL, whether_abort_operation]];   

Pos.whether_abort_operation = whether_abort_operation;

if ~whether_abort_operation
    % Copying data to output data structure
    Pos.trueRange = range;
    Pos.rangeResid = Res;
    Pos.nrSats = diff([0 nrSatsUsed]);
    Pos.signals = allSettings.sys.enabledSignals;
    Pos.xyz = pos(1:3);
    Pos.dt = pos(4:end);
    Pos.W = W;
    
    % % Get dop values
    % Pos.dop = getDOPValues(allSettings.const,H, Pos.xyz);
    % 
    % Calculate fom
    Pos.fom = norm(Res/length(Res));
    
    % Check if solution is valid
    if(Pos.fom < 50)
        Pos.bValid = true;
    end
else
    % Calculate fom (set a large value)
    Pos.fom = 1e5;
end

if allSettings.sys.raim_type ~= 2 % having RAIM
    if whether_compute_PL
        PL = compute_PL(Updated_mat, Detect_results, allSettings);
    else
        PL = NaN;
    end
    Pos.raim = [Detect_results.WSSE_sqrt, Detect_results.Thres, PL, whether_abort_operation];
end

 
 


