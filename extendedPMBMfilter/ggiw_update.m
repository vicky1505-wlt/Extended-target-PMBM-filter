function [ Log_L,alpha_upd,beta_upd,x_upd,P_upd,v_upd,V_upd ] = ggiw_update( H,meas,alpha,beta,x,P,v,V,model )
% GGIW update

d = 2;
R = model.R;
W = size(meas,2);

% update Gamma distribution parameters
alpha_upd = alpha + W;
beta_upd = beta + 1;

X_hat = V./(v-2*d-2);
R_hat = X_hat + R;
sqrtRhat = sqrtm(R_hat);
S = H*P*H' + R_hat./W;
K = P*H'/S;

% update Gaussian distribution parameters
z_bar = mean(meas,2);
epsilon = z_bar - H*x;
x_upd = x + K*epsilon;
P_upd = P - K*H*P;

% update inverse Wishart distribution parameters
Z = (meas - z_bar)*(meas - z_bar)';
N = epsilon*epsilon';
Xsqrt = sqrtm(X_hat);
Ssqrt = sqrtm(S);
N_hat = (Xsqrt/Ssqrt)*N*(Ssqrt\Xsqrt);
Zhat = (Xsqrt/sqrtRhat)*Z*(sqrtRhat\Xsqrt);
V_upd = V + N_hat + Zhat;
v_upd = v + W;

% Compute GGIW predictive likelihood

temp1 = (v-d-1)/2;
temp2 = (v_upd-d-1)/2;

log_temp3 = temp1*log(det(V))-temp2*log(det(V_upd));
log_temp4 = logamma2(temp2)-logamma2(temp1);
log_temp5 = 1/2*log(det(X_hat))-1/2*log(det(S))-1/2*log(det(R_hat));
log_temp6 = gammaln(alpha_upd)-gammaln(alpha);
log_temp7 = alpha*log(beta)-alpha_upd*log(beta_upd);

Log_L = (-d/2)*(W*log(pi)+log(W))+log_temp3+log_temp4+log_temp5+log_temp6+log_temp7;

end

