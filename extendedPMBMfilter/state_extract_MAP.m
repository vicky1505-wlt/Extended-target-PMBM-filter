function [ est ] = state_extract_MAP( ggiw_mb )

d = 2;
r = ggiw_mb.r;
r(r>1-eps) = 1-eps; % minus a small number to ensure that 1-r>0
ss = false(size(r));
pcard = prod(1-r)*poly(-r./(1-r));
[~,n] = max(pcard);
[~,o] = sort(-r);
n = n - 1;
ss(o(1:n)) = true;
est.x = ggiw_mb.x(:,ss);
v = ggiw_mb.v(ss);
V = ggiw_mb.V(:,:,ss);

for i = 1:length(v)
    est.X(:,:,i) = V(:,:,i)./(v(i)-2*d-2);
end
if isempty(v)
    est.X = zeros(2,2,0);
end

end