function [ggiw_mb,ggiw_ppp] = predicting(ggiw_mb,ggiw_ppp,model)

F = model.F;
Q = model.Q;
Ps = model.Ps;
Ts = model.Ts;
eta = model.eta;
tao = model.tao;
d = 2;
M = eye(2);

alpha_b = model.alpha_b;
beta_b = model.beta_b;
xb = model.xb;
Pb = model.Pb;
vb = model.vb;
Vb = model.Vb;
wb = model.wb;

% GGIW-MBM predict
J = length(ggiw_mb);
for j = 1:J
    ggiw_mb{j}.r = Ps*ggiw_mb{j}.r;
    ggiw_mb{j}.alpha = ggiw_mb{j}.alpha/eta;
    ggiw_mb{j}.beta = ggiw_mb{j}.beta/eta;
    ggiw_mb{j}.x = F*ggiw_mb{j}.x;
    n = length(ggiw_mb{j}.r);
    for i = 1:n
        ggiw_mb{j}.P(:,:,i) = F*ggiw_mb{j}.P(:,:,i)*F' + Q;
        ggiw_mb{j}.v(i) = 2*d+2+exp(-Ts/tao)*(ggiw_mb{j}.v(i)-2*d-2);
        ggiw_mb{j}.V(:,:,i) = exp(-Ts/tao)*M*ggiw_mb{j}.V(:,:,i)*M';
    end
end

% GGIW-PPP predict
ggiw_ppp.wu = [wb;ggiw_ppp.wu*Ps];
ggiw_ppp.alpha_u = ggiw_ppp.alpha_u/eta;
ggiw_ppp.beta_u = ggiw_ppp.beta_u/eta;
ggiw_ppp.xu = F*ggiw_ppp.xu;
n = length(ggiw_ppp.alpha_u);
for i = 1:n
    ggiw_ppp.Pu(:,:,i) = F*ggiw_ppp.Pu(:,:,i)*F' + Q;
    ggiw_ppp.vu(i) = 2*d+2+exp(-Ts/tao)*(ggiw_ppp.vu(i)-2*d-2);
    ggiw_ppp.Vu(:,:,i) = exp(-Ts/tao)*M*ggiw_ppp.Vu(:,:,i)*M';
end
ggiw_ppp.alpha_u = [alpha_b;ggiw_ppp.alpha_u];
ggiw_ppp.beta_u = [beta_b;ggiw_ppp.beta_u];
ggiw_ppp.xu = [xb ggiw_ppp.xu];
ggiw_ppp.Pu = cat(3,Pb,ggiw_ppp.Pu);
ggiw_ppp.vu = [vb;ggiw_ppp.vu];
ggiw_ppp.Vu = cat(3,Vb,ggiw_ppp.Vu);

end

