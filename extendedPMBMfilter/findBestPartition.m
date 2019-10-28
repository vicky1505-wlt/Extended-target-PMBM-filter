function [ bestPartition ] = findBestPartition( ggiw_ppp,Partitions,z,model )
% Find the best partition

if isempty(Partitions)
    bestPartition = cell(0,1);
    return;
end

H = model.H;
Pd = model.Pd;
lambda_fa = model.lambda_fa;

P = length(Partitions);         % number of partitions
nGGIW = length(ggiw_ppp.wu);    % number of GGIWs in PPP

log_Lu = zeros(nGGIW,1);
lp = zeros(P,1);
for p = 1:P
    np = length(Partitions{p}); % number of cells in each partition
    lc = zeros(np,1);
    
    for i = 1:np
        c = z(:,Partitions{p}{i});
        for u = 1:nGGIW
            % compute predictive likelihood
            log_Lu(u) = ggiw_update(H,c,ggiw_ppp.alpha_u(u),ggiw_ppp.beta_u(u),...
                ggiw_ppp.xu(:,u),ggiw_ppp.Pu(:,:,u),ggiw_ppp.vu(u),ggiw_ppp.Vu(:,:,u),model);
        end
        lc(i) = Pd*ggiw_ppp.wu'*exp(log_Lu);
        % for measurement cell single measurement, also considers the
        % likelihood of clutter
        if length(Partitions{p}{i}) == 1
            lc(i) = lc(i) + lambda_fa;
        end
    end
    % joint likelihood for a measurement partition
    lp(p) = prod(lc);
end

% select the most likely partition
[~,p_best] = max(lp);
bestPartition = Partitions{p_best};

end

