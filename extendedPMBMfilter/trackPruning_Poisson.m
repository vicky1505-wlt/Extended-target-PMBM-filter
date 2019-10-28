function [ ggiw_ppp,idx ] = trackPruning_Poisson( ggiw_ppp, model )

d = 2;
idx1 = ggiw_ppp.wu > model.threshold_u;
idx2 = -log((ggiw_ppp.beta_u./(ggiw_ppp.beta_u+1)).^ggiw_ppp.alpha_u)>1;
idx3 = ggiw_ppp.vu > 2*d + 2;

idx = idx1&idx2&idx3;

ggiw_ppp.wu = ggiw_ppp.wu(idx);
ggiw_ppp.alpha_u = ggiw_ppp.alpha_u(idx);
ggiw_ppp.beta_u = ggiw_ppp.beta_u(idx);
ggiw_ppp.xu = ggiw_ppp.xu(:,idx);
ggiw_ppp.Pu = ggiw_ppp.Pu(:,:,idx);
ggiw_ppp.vu = ggiw_ppp.vu(idx);
ggiw_ppp.Vu = ggiw_ppp.Vu(:,:,idx);

end

