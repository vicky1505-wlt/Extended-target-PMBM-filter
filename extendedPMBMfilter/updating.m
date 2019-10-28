function [ggiw_mbm_update,ggiw_ppp_update,W,est] = updating(ggiw_mbm,ggiw_ppp,Z,model,W)

ggiw_mbm_update = cell(0,1);
J = length(ggiw_mbm);
w = cell(J,1);

% PPP update
ggiw_ppp_update = ppp_update(ggiw_ppp,model);

if isempty(ggiw_mbm{1}.r)
    Z_gate = gate_meas_gms(Z,ggiw_ppp,model);
    partitions = genPartitions_dbscan(Z_gate,model);
    bestPartition = findBestPartition(ggiw_ppp,partitions,Z_gate,model);
    [ggiw_new,~] = newBernoulli(ggiw_ppp,bestPartition,Z_gate,model);
    ggiw_mbm_update{1} = trackPruning_Bernoulli(ggiw_new,model);
    est = state_extract_MBM(W,ggiw_mbm_update);
    return;
end

% Loop through each global hypothesis
for j = 1:J
    
    % Gating and clustering for pre-existing targets
    [gatingGroup,gatingGGIW,idx_out] = groupGating(Z,ggiw_mbm{j},model);
    
    % Gating for unknown targets; find the measurements that are outside
    % the gates of pre-existing targets but are inside the gates of unknown
    % targets
    Z_gate = gate_meas_gms(Z(:,idx_out),ggiw_ppp,model);
    
    % Generate measurement partitions using db-scan
    partitions = genPartitions_dbscan(Z_gate,model);
    
    % Find the most likely measurement partition for unknown targets
    bestPartition = findBestPartition(ggiw_ppp,partitions,Z_gate,model);
    
    % Create the new GGIW components for unknown targets
    [ggiw_new,Lc] = newBernoulli(ggiw_ppp,bestPartition,Z_gate,model);
    
    % Obtain Mt-best data association hypotheses using Murty's algorithm
    Mt = ceil(model.M*W(j));
    
    % Perform MBM update for cluster independently
    G = length(gatingGroup);
    ggiw_mbm_gate = cell(G,1);
    Wg = cell(G,1);
    
    for g = 1:G
        
        % Generate measurement partitions using db-scan
        MPs = genPartitions_dbscan(Z(:,gatingGroup{g}),model);
        nP = length(MPs);   
        if isempty(MPs) % this corresponds to a cluster w/o measurement
            nP = 1;
            MPs = cell(1,1);
        end
        
        log_Wp = zeros(0,1);
        ggiw_mbm_gate_hat = cell(0,1);
        for i = 1:nP
            % Create Bernoulli components for PPP
            [ggiw_mb_new,Lnew] = newBernoulli(ggiw_ppp,MPs{i},Z(:,gatingGroup{g}),model);
            
            % Create Bernoulli components for MB
            [ggiw_mb_miss,ggiw_mb_upd,Lmiss,Lupd] = mbm_update(gatingGGIW{g},model,Z(:,gatingGroup{g}),MPs{i});
            
            % Data association
            [ggiw_mbm_murty,log_Wmurty] = dataAssociation(ggiw_mb_new,Lnew,...
                ggiw_mb_upd,Lupd,ggiw_mb_miss,Lmiss,model,Mt);
            
            % Concatenate MBM obtained for each MB
            log_Wp = [log_Wp;log_Wmurty];
            ggiw_mbm_gate_hat = cat(1,ggiw_mbm_gate_hat,ggiw_mbm_murty);
        end
        
        % Pruning by only keeping MBs with total weights that correspond to
        % 1-model.threshold_w
        Wp_normalized = exp(normalizeLogWeights(log_Wp));
        [Wp_sorted,order] = sort(Wp_normalized,'descend');
        pos = find(cumsum(Wp_sorted) >= 1-model.threshold_w,1);
        Wg{g} = log_Wp(order(1:pos));
        ggiw_mbm_gate{g} = ggiw_mbm_gate_hat(order(1:pos));

    end
    
    % Combine MBM from different gating groups
    [ggiw_mbm_upd,Wupd] = hypo_update(Wg,ggiw_mbm_gate,ggiw_new,model);
    
    % Concatenate MBs for different global hypotheses
    ggiw_mbm_update = cat(1,ggiw_mbm_update,ggiw_mbm_upd);
    w{j} = log(W(j)) + Wupd + sum(log(Lc));
end

% Normalise global hypotheses weights
W = cell2mat(w);
W = exp(normalizeLogWeights(W));

% Prune low-weight global hypothesis
[W_sorted,order] = sort(W,'descend');
pos = find(cumsum(W_sorted) >= 1-model.threshold_w,1);
W = W(order(1:pos));
ggiw_mbm_update = ggiw_mbm_update(order(1:pos));
W = W/sum(W);

% MBM Merging, merge similar MBs in the sense of small KL-divergence
[ggiw_mbm_update,log_Wmbm_m] = MergeMBM(ggiw_mbm_update,log(W),model);
W = exp(normalizeLogWeights(log_Wmbm_m));

% Recycling
[ggiw_mbm_update,ggiw_ppp_update] = recycling_ett(ggiw_mbm_update,ggiw_ppp_update,model);

% Target states extraction
est = state_extract_MBM(W,ggiw_mbm_update);

% In case, all the MBs are pruned; if happens, initialise the parameter
% structure again
if isempty(ggiw_mbm_update)
    W = 1;
    % Initial GGIW-MB distribution
    x_dim = model.x_dim;
    z_dim = model.z_dim;
    % each cell corresponds to a GGIW-MB
    ggiw_mbm_update{1}.r = zeros(0,1);
    ggiw_mbm_update{1}.x = zeros(x_dim,0);
    ggiw_mbm_update{1}.P = zeros(x_dim,x_dim,0);
    ggiw_mbm_update{1}.alpha = zeros(0,1);
    ggiw_mbm_update{1}.beta = zeros(0,1);
    ggiw_mbm_update{1}.v = zeros(0,1);
    ggiw_mbm_update{1}.V = zeros(z_dim,z_dim,0);
end

end

