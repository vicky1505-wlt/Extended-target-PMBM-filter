function [ ggiw_mbm_murty,log_nW ] = dataAssociation( ggiw_mb_new,...
    Lnew,ggiw_mb_upd,Lupd,ggiw_mb_miss,Lmiss,model,Mt)

dim = model.x_dim;
n = size(Lupd,1);
m = size(Lupd,2);
cost = -log(Lupd./repmat(Lnew',n,1));

% Obtain M-best data association hypotheses using Murty's algorithm
[bestAssign,nCost] = mbestwrap_updt_custom(cost,Mt,Lmiss);

log_nW = -nCost'+sum(log(Lnew));

M = length(nCost);
ggiw_mbm_murty = cell(M,1);
for j = 1:M
    
    r = zeros(n+m,1);
    alpha = zeros(n+m,1);
    beta = zeros(n+m,1);
    x = zeros(dim,n+m);
    P = zeros(dim,dim,n+m);
    v = zeros(n+m,1);
    V = zeros(2,2,n+m);
    for i=1:n
        if bestAssign(j,i) ~= 0
            % measurement update for pre-existing targets
            r(i) = ggiw_mb_upd.r(i,bestAssign(j,i));
            alpha(i) = ggiw_mb_upd.alpha(i,bestAssign(j,i));
            beta(i) = ggiw_mb_upd.beta(i,bestAssign(j,i));
            x(:,i) = ggiw_mb_upd.x(:,i,bestAssign(j,i));
            P(:,:,i) = ggiw_mb_upd.P(:,:,i,bestAssign(j,i));
            v(i) = ggiw_mb_upd.v(i,bestAssign(j,i));
            V(:,:,i) = ggiw_mb_upd.V(:,:,i,bestAssign(j,i));
        else
            % misdetection update
            r(i) = ggiw_mb_miss.r(i);
            alpha(i) = ggiw_mb_miss.alpha(i);
            beta(i) = ggiw_mb_miss.beta(i);
            x(:,i) = ggiw_mb_miss.x(:,i);
            P(:,:,i) = ggiw_mb_miss.P(:,:,i);
            v(i) = ggiw_mb_miss.v(i);
            V(:,:,i) = ggiw_mb_miss.V(:,:,i);
        end
    end
    
    for i=n+1:n+m
        if ~any(bestAssign(j,:)==i-n)
            % measurement update for unknown targets
            r(i) = ggiw_mb_new.r(i-n);
            alpha(i) = ggiw_mb_new.alpha(i-n);
            beta(i) = ggiw_mb_new.beta(i-n);
            x(:,i) = ggiw_mb_new.x(:,i-n);
            P(:,:,i) = ggiw_mb_new.P(:,:,i-n);
            v(i) = ggiw_mb_new.v(i-n);
            V(:,:,i) = ggiw_mb_new.V(:,:,i-n);
        end
    end 
    
    ggiw_mb.r = r;
    ggiw_mb.alpha = alpha;
    ggiw_mb.beta = beta;
    ggiw_mb.x = x;
    ggiw_mb.P = P;
    ggiw_mb.v = v;
    ggiw_mb.V = V;
    
    ggiw_mbm_murty{j} = ggiw_mb;
    
end


end

