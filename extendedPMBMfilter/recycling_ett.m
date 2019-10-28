function [ ggiw_mbm_update,ggiw_ppp_update ] = recycling_ett( ggiw_mbm_update,ggiw_ppp_update,model)

% Recycling: approximate Bernoulli components with very small existence
% probability as PPP

G = length(ggiw_mbm_update);
ggiw_mbm_loweight.r = zeros(0,1);
ggiw_mbm_loweight.alpha = zeros(0,1);
ggiw_mbm_loweight.beta = zeros(0,1);
ggiw_mbm_loweight.x = zeros(4,0);
ggiw_mbm_loweight.P = zeros(4,4,0);
ggiw_mbm_loweight.v = zeros(0,1);
ggiw_mbm_loweight.V = zeros(2,2,0);

for g = 1:G
    idx = ggiw_mbm_update{g}.r >= model.threshold_recycle;
    ggiw_mbm_loweight.r = [ggiw_mbm_loweight.r;ggiw_mbm_update{g}.r(~idx)];
    ggiw_mbm_loweight.alpha = [ggiw_mbm_loweight.alpha;ggiw_mbm_update{g}.alpha(~idx)];
    ggiw_mbm_loweight.beta = [ggiw_mbm_loweight.beta;ggiw_mbm_update{g}.beta(~idx)];
    ggiw_mbm_loweight.x = [ggiw_mbm_loweight.x ggiw_mbm_update{g}.x(:,~idx)];
    ggiw_mbm_loweight.P = cat(3,ggiw_mbm_loweight.P,ggiw_mbm_update{g}.P(:,:,~idx));
    ggiw_mbm_loweight.v = [ggiw_mbm_loweight.v;ggiw_mbm_update{g}.v(~idx)];
    ggiw_mbm_loweight.V = cat(3,ggiw_mbm_loweight.V,ggiw_mbm_update{g}.V(:,:,~idx));
    ggiw_mbm_update{g}.r = ggiw_mbm_update{g}.r(idx);
    ggiw_mbm_update{g}.alpha = ggiw_mbm_update{g}.alpha(idx);
    ggiw_mbm_update{g}.beta = ggiw_mbm_update{g}.beta(idx);
    ggiw_mbm_update{g}.x = ggiw_mbm_update{g}.x(:,idx);
    ggiw_mbm_update{g}.P = ggiw_mbm_update{g}.P(:,:,idx);
    ggiw_mbm_update{g}.v = ggiw_mbm_update{g}.v(idx);
    ggiw_mbm_update{g}.V = ggiw_mbm_update{g}.V(:,:,idx);
end

% Truncation
l = 0;
I = find(ggiw_mbm_loweight.r > model.threshold_r & ...
    -log((ggiw_mbm_loweight.beta./(ggiw_mbm_loweight.beta+1)).^ggiw_mbm_loweight.alpha) > 1 & ...
    ggiw_mbm_loweight.v > 6);
if isempty(I)
    ggiw_mb_hat.r = zeros(0,1);
    ggiw_mb_hat.alpha = zeros(0,1);
    ggiw_mb_hat.beta = zeros(0,1);
    ggiw_mb_hat.x = zeros(4,0);
    ggiw_mb_hat.P = zeros(4,4,0);
    ggiw_mb_hat.v = zeros(0,1);
    ggiw_mb_hat.V = zeros(2,2,0);
end

while (~isempty(I))
    l = l + 1;
    [~,idx] = max(ggiw_mbm_loweight.r(I));
    L = [];
    for i = 1:length(I)
        temp = GGIW_KLdiff2(ggiw_mbm_loweight.alpha(I(idx)),...
            ggiw_mbm_loweight.beta(I(idx)),ggiw_mbm_loweight.x(:,I(idx)),...
            ggiw_mbm_loweight.P(:,:,I(idx)),ggiw_mbm_loweight.v(I(idx)),...
            ggiw_mbm_loweight.V(:,:,I(idx)),ggiw_mbm_loweight.alpha(I(i)),...
            ggiw_mbm_loweight.beta(I(i)),ggiw_mbm_loweight.x(:,I(i)),...
            ggiw_mbm_loweight.P(:,:,I(i)),ggiw_mbm_loweight.v(I(i)),ggiw_mbm_loweight.V(:,:,I(i)));
        if temp <= model.mergingThreshold_recycle
            L = [L I(i)];
        end
    end
    [ggiw_mb_hat.r(l),ggiw_mb_hat.alpha(l),ggiw_mb_hat.beta(l),ggiw_mb_hat.x(:,l),ggiw_mb_hat.P(:,:,l),...
        ggiw_mb_hat.v(l),ggiw_mb_hat.V(:,:,l)] = GGIW_merge(ggiw_mbm_loweight.r(L),...
        ggiw_mbm_loweight.alpha(L),ggiw_mbm_loweight.beta(L),ggiw_mbm_loweight.x(:,L),...
        ggiw_mbm_loweight.P(:,:,L),ggiw_mbm_loweight.v(L),ggiw_mbm_loweight.V(:,:,L));
    I = setdiff(I,L);
end

ggiw_ppp_update = recycle(ggiw_ppp_update,ggiw_mb_hat);

end

function ggiw_ppp_update = recycle(ggiw_ppp_update,ggiw_mb_hat)
ggiw_ppp_update.wu = [ggiw_ppp_update.wu;ggiw_mb_hat.r'];
ggiw_ppp_update.alpha_u = [ggiw_ppp_update.alpha_u;ggiw_mb_hat.alpha'];
ggiw_ppp_update.beta_u = [ggiw_ppp_update.beta_u;ggiw_mb_hat.beta'];
ggiw_ppp_update.xu = [ggiw_ppp_update.xu ggiw_mb_hat.x];
ggiw_ppp_update.Pu = cat(3,ggiw_ppp_update.Pu,ggiw_mb_hat.P);
ggiw_ppp_update.vu = [ggiw_ppp_update.vu;ggiw_mb_hat.v'];
ggiw_ppp_update.Vu = cat(3,ggiw_ppp_update.Vu,ggiw_mb_hat.V);
end

