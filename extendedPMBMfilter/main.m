dbstop if error
clc;clear

load ScenarioWith27targets.mat

para.false_alarm_rate = 60; % choose false alarm rate
para.detection_prob = 0.9;  % choose detection probability
 
% Model paramters and ground truth
[groundTruth,Z,model] = gen_data_many_targets(para,targetTracks);

% Total time step
K = length(groundTruth);

% (G)OSPA parameters
gospa_p = 1;
gospa_c = 10;
gospa_vals = zeros(K,4);
ospa_vals = zeros(K,3);


estimationResults = cell(K,1);

% Initialise multi-Bernoulli mixture and PPP parameters
[ggiw_mbm,ggiw_ppp] = Initialisation(model);
W = 1;

fprintf('Current time step: ')
for t = 1:K
    
    fprintf('%g ', t);
    
    % Prediction step
    [ggiw_mbm,ggiw_ppp] = predicting(ggiw_mbm,ggiw_ppp,model);
    
    % Update step
    [ggiw_mbm,ggiw_ppp,W,est] = updating(ggiw_mbm,ggiw_ppp,Z{t},model,W);
    
    % Performance evaluation using GOSPA
    gospa_vals(t,:) = GOSPAmetric(est,groundTruth{t},gospa_c,gospa_p);
    % Performance evaluation using OSPA
    ospa_vals(t,:) = OSPAmetric(est,groundTruth{t},gospa_c,gospa_p);
    
    % Store estimates
    estimationResults{t} = est;
end

% Compute (G)OSPA error, averaged per MB run and per time step
averGospa = mean(mean(gospa_vals,3));
averOspa = mean(mean(ospa_vals,3));

fprintf('\nGOSPA error: %g, Localisation error: %g, Missed detection: %g, False detection: %g', ...
    averGospa(1), averGospa(2), averGospa(3), averGospa(4));

fprintf('\nOSPA error: %g, Localisation error: %g, Cardinality error: %g\n', ...
    averOspa(1), averOspa(2), averOspa(3));
