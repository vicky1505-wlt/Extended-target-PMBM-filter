function [ ggiw_mbm, ggiw_ppp ] = Initialisation( model )

% Initial GGIW-MB distribution
x_dim = model.x_dim;
z_dim = model.z_dim;
% each cell corresponds to a GGIW-MB
ggiw_mbm{1}.r = zeros(0,1);
ggiw_mbm{1}.x = zeros(x_dim,0);
ggiw_mbm{1}.P = zeros(x_dim,x_dim,0);
ggiw_mbm{1}.alpha = zeros(0,1);
ggiw_mbm{1}.beta = zeros(0,1);
ggiw_mbm{1}.v = zeros(0,1);
ggiw_mbm{1}.V = zeros(z_dim,z_dim,0);

% Initial unknown target GGIW-PPP parameters
ggiw_ppp.alpha_u = model.alpha_b;
ggiw_ppp.beta_u = model.beta_b;
ggiw_ppp.xu = model.xb;
ggiw_ppp.Pu = model.Pb;
ggiw_ppp.vu = model.vb;
ggiw_ppp.Vu = model.Vb;
ggiw_ppp.wu = model.wb;

end

