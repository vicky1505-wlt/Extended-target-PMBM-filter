function [ ggiw_mb,idx ] = trackPruning_Bernoulli( ggiw_mb, model )

d = 2;
idx1 = ggiw_mb.r > model.threshold_r;
idx2 = -log((ggiw_mb.beta./(ggiw_mb.beta+1)).^ggiw_mb.alpha)>1;
idx3 = ggiw_mb.v > 2*d + 2;

idx = idx1&idx2&idx3;

ggiw_mb.r = ggiw_mb.r(idx);
ggiw_mb.alpha = ggiw_mb.alpha(idx);
ggiw_mb.beta = ggiw_mb.beta(idx);
ggiw_mb.x = ggiw_mb.x(:,idx);
ggiw_mb.P = ggiw_mb.P(:,:,idx);
ggiw_mb.v = ggiw_mb.v(idx);
ggiw_mb.V = ggiw_mb.V(:,:,idx);

end

