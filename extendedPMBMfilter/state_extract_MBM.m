function [ est ] = state_extract_MBM( W, ggiw_mbm )

if isempty(W)
    est.x = zeros(4,0);
    est.X = zeros(2,2,0);
    return;
end

[~,idx] = max(W);
ggiw_mb = ggiw_mbm{idx};

d = 2;
% select Bernoulli density with existence probability larger than 0.5
ss = ggiw_mb.r >= 0.5;
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