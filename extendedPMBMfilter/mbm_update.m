function [ggiw_mb_miss,ggiw_mb_upd,Lmiss,Lupd] = mbm_update(ggiw_mb,model,Z,Partition)

d = 2;
Pd = model.Pd;
x_dim = model.x_dim;
H = model.H;
np = length(Partition);         % number of measurement cells

nBernoulli = length(ggiw_mb.r); % number of Bernoullis
temp = (ggiw_mb.beta./(ggiw_mb.beta+1)).^ggiw_mb.alpha;
qD = 1 - Pd + Pd*temp;          % missed detection likelihood

% Missed detection hypothesis
Lmiss = 1 - ggiw_mb.r + ggiw_mb.r.*qD;
ggiw_mb_miss.r = ggiw_mb.r.*qD./Lmiss;
w_miss1 = (1-Pd)./qD;
w_miss2 = Pd*temp./qD;

% Merge two GGIW components created for missed detection hypotheses
ggiw_mb_miss.alpha = ggiw_mb.alpha;
beta = ggiw_mb.beta;
ggiw_mb_miss.beta = zeros(nBernoulli,1);
for i = 1:nBernoulli
    w_bar = w_miss1(i) + w_miss2(i);
    temp = 1/w_bar*ggiw_mb_miss.alpha(i)*(w_miss1(i)/beta(i)+w_miss2(i)/(beta(i)+1));
    beta_hat = ggiw_mb_miss.alpha(i)/temp;
    ggiw_mb_miss.beta(i) = beta_hat;
end
ggiw_mb_miss.x = ggiw_mb.x;
ggiw_mb_miss.P = ggiw_mb.P;
ggiw_mb_miss.v = ggiw_mb.v;
ggiw_mb_miss.V = ggiw_mb.V;

% Measurement-updated hypothesis
ggiw_mb_upd.r = ones(nBernoulli,np);
ggiw_mb_upd.alpha = zeros(nBernoulli,np);
ggiw_mb_upd.beta = zeros(nBernoulli,np);
ggiw_mb_upd.x = zeros(x_dim,nBernoulli,np);
ggiw_mb_upd.P = zeros(x_dim,x_dim,nBernoulli,np);
ggiw_mb_upd.v = zeros(nBernoulli,np);
ggiw_mb_upd.V = zeros(d,d,nBernoulli,np);
Lupd = zeros(nBernoulli,np);

for i = 1:nBernoulli
    for j = 1:np
        [log_L,ggiw_mb_upd.alpha(i,j),ggiw_mb_upd.beta(i,j),ggiw_mb_upd.x(:,i,j),...
            ggiw_mb_upd.P(:,:,i,j),ggiw_mb_upd.v(i,j),ggiw_mb_upd.V(:,:,i,j)] = ...
            ggiw_update(H,Z(:,Partition{j}),ggiw_mb.alpha(i),ggiw_mb.beta(i),ggiw_mb.x(:,i),...
            ggiw_mb.P(:,:,i),ggiw_mb.v(i),ggiw_mb.V(:,:,i),model);
        % Data association likelihood
        Lupd(i,j) = ggiw_mb.r(i)*exp(log_L)*model.Pd;
    end
end

end

