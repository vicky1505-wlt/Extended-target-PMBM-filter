function [ ggiw_new_merge,Lc ] = newBernoulli( ggiw_ppp,Partition,z,model )
% GGIW PPP UPDATE: DETECTION

H = model.H;
Pd = model.Pd;
lambda_fa = model.lambda_fa;
x_dim = model.x_dim;
d = 2;

% new Bernoulli estimation for the selected partition
np = length(Partition);       % number of cells
nGGIW = length(ggiw_ppp.wu);  % number of GGIWs
log_Lu = -inf(nGGIW,np);
Lc = zeros(np,1);
rc = ones(np,1);
ggiw_new = cell(np,1);

for i = 1:np
    
    W = length(Partition{i}); % number of measurements in each cell
    c = z(:,Partition{i});
    ggiw_new{i}.alpha_u = zeros(nGGIW,1);
    ggiw_new{i}.beta_u = zeros(nGGIW,1);
    ggiw_new{i}.vu = zeros(nGGIW,1);
    ggiw_new{i}.xu = zeros(4,nGGIW);
    ggiw_new{i}.Pu = zeros(4,4,nGGIW);
    ggiw_new{i}.Vu = zeros(2,2,nGGIW);
    ggiw_new{i}.wu = zeros(nGGIW,1);
    
    if ~isempty(gate_meas_gms(c,ggiw_ppp,model))
        for u = 1:nGGIW
            [log_Lu(u,i),ggiw_new{i}.alpha_u(u),ggiw_new{i}.beta_u(u),ggiw_new{i}.xu(:,u),...
                ggiw_new{i}.Pu(:,:,u),ggiw_new{i}.vu(u),ggiw_new{i}.Vu(:,:,u)] = ...
                ggiw_update(H,c,ggiw_ppp.alpha_u(u),ggiw_ppp.beta_u(u),...
                ggiw_ppp.xu(:,u),ggiw_ppp.Pu(:,:,u),ggiw_ppp.vu(u),ggiw_ppp.Vu(:,:,u),model);
        end
        log_Lu(:,i) = log_Lu(:,i) + log(Pd);
        Lc(i) = ggiw_ppp.wu'*exp(log_Lu(:,i));
        if W == 1
            % for measurement cell with single measurement, also considers
            % the likelihood of clutter
            rc(i) = Lc(i)/(Lc(i)+lambda_fa);
            log_Lu(:,i) = log_Lu(:,i) + log(lambda_fa);
        end
        Lc(i) = ggiw_ppp.wu'*exp(log_Lu(:,i));
        ggiw_new{i}.wu = ggiw_ppp.wu.*exp(log_Lu(:,i))/Lc(i);
    end
    
end

ggiw_new_merge.r = rc;
ggiw_new_merge.alpha = zeros(np,1);
ggiw_new_merge.beta = zeros(np,1);
ggiw_new_merge.x = zeros(x_dim,np);
ggiw_new_merge.P = zeros(x_dim,x_dim,np);
ggiw_new_merge.v = zeros(np,1);
ggiw_new_merge.V = zeros(d,d,np);
% GGIW mixture reduction
for i = 1:np
    valid_idx = ggiw_new{i}.wu > 1e-3;
    if any(valid_idx)
        [~,ggiw_new_merge.alpha(i),ggiw_new_merge.beta(i),ggiw_new_merge.x(:,i),...
            ggiw_new_merge.P(:,:,i),ggiw_new_merge.v(i),ggiw_new_merge.V(:,:,i)] = ...
            GGIW_merge(ggiw_new{i}.wu(valid_idx),ggiw_new{i}.alpha_u(valid_idx),ggiw_new{i}.beta_u(valid_idx),...
            ggiw_new{i}.xu(:,valid_idx),ggiw_new{i}.Pu(:,:,valid_idx),ggiw_new{i}.vu(valid_idx),ggiw_new{i}.Vu(:,:,valid_idx));
    else
        ggiw_new_merge.r(i) = 0;
        Lc(i) = eps;
    end
end

end