function [z_gate]= gate_meas_gms(z,ggiw_ppp,model)

d = 2;
H = model.H;
R = model.R;
gamma = model.gamma;
x = ggiw_ppp.xu;
P = ggiw_ppp.Pu;
v = ggiw_ppp.vu;
V = ggiw_ppp.Vu;

zlength = size(z,2);
plength = size(x,2);

in_gate = false(zlength,1);

for j=1:plength
    Sj= V(:,:,j)/(v(j)-2*d-2) + H*P(:,:,j)*H' + R;
    nu= z- model.H*repmat(x(:,j),[1 zlength]);
    dist= sum((inv(chol(Sj))'*nu).^2);
    in_gate(dist<gamma) = true;
end
z_gate = z(:,in_gate);