function [ ggiw_mbm_upd, Wupd ] = hypo_update( Wg, ggiw_mbm_gate, ggiw_new, model )

% Pruning unlikely single target hypotheses
ggiw_new = trackPruning_Bernoulli(ggiw_new,model);

N = length(Wg);
for i = 1:N
    for g = 1:length(ggiw_mbm_gate{i})
        % prune unlikely single targets hypotheses
        ggiw_mbm_gate{i}{g} = trackPruning_Bernoulli(ggiw_mbm_gate{i}{g},model);
    end
end

% get all the possible different combinations of MBM
sets = cell(1,N);
for i = 1:N
    sets{i} = 1:length(Wg{i});
end

if ~isempty(sets)
    c = allcombs(sets{:});
    [N,M] = size(c);
else
    N = 0;
    M = 0;
end

ggiw_mbm_upd = cell(N,1);
Wupd = zeros(N,1);
for i = 1:N
    r = zeros(1,0);
    alpha = zeros(1,0);
    beta = zeros(1,0);
    x = zeros(4,0);
    P = zeros(4,4,0);
    v = zeros(1,0);
    V = zeros(2,2,0);
    wc = 0;
    
    for k = 1:M
        wc = wc+Wg{k}(c(i,k));
        r = [r;ggiw_mbm_gate{k}{c(i,k)}.r];
        alpha = [alpha;ggiw_mbm_gate{k}{c(i,k)}.alpha];
        beta = [beta;ggiw_mbm_gate{k}{c(i,k)}.beta];
        x = [x ggiw_mbm_gate{k}{c(i,k)}.x];
        P = cat(3,P,ggiw_mbm_gate{k}{c(i,k)}.P);
        v = [v;ggiw_mbm_gate{k}{c(i,k)}.v];
        V = cat(3,V,ggiw_mbm_gate{k}{c(i,k)}.V);
    end
    Wupd(i) = wc;
    
    ggiw_mbm_upd{i}.r = [r;ggiw_new.r];
    ggiw_mbm_upd{i}.alpha = [alpha;ggiw_new.alpha];
    ggiw_mbm_upd{i}.beta = [beta;ggiw_new.beta];
    ggiw_mbm_upd{i}.x = [x ggiw_new.x];
    ggiw_mbm_upd{i}.P = cat(3,P,ggiw_new.P);
    ggiw_mbm_upd{i}.v = [v;ggiw_new.v];
    ggiw_mbm_upd{i}.V = cat(3,V,ggiw_new.V);
    
end

end

