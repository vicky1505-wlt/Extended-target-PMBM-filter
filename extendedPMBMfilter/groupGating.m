function [gatingGroup,gatingGGIW,idx_out] = groupGating(z,ggiw_mb,model)

d = 2;
x_dim = model.x_dim;
H = model.H;
R = model.R;
gamma = model.gamma;
x = ggiw_mb.x;
P = ggiw_mb.P;
v = ggiw_mb.v;
V = ggiw_mb.V;

idx_in = [];
zlength = size(z,2);
plength = size(x,2);
tentativeGroup = cell(plength,1);

% Form tentative gating groups
for j=1:plength
    Sj = V(:,:,j)/(v(j)-2*d-2) + H*P(:,:,j)*H' + R;
    Vs= chol(Sj);
    inv_sqrt_Sj= inv(Vs);
    nu = z - H*repmat(x(:,j),[1 zlength]);
    dist = sum((inv_sqrt_Sj'*nu).^2);
    idx_valid = find(dist < gamma);
    idx_in= union(idx_in,idx_valid);
    tentativeGroup{j} = idx_valid;
end

% Get measurement indices out of the gate
idx_out = setdiff(1:zlength,idx_in);

% Merging gating groups with common measurements
label = 1:plength;
for i = 1:plength-1
    for j = i+1:plength
        if ~isempty(intersect(tentativeGroup{i},tentativeGroup{j}))
            label(j) = min(label(i),label(j));
        end
    end
end

% Clustering
uniqueLabel = unique(label);
n = length(uniqueLabel);
gatingGroup = cell(n,1);
gatingGGIW = cell(n,1);
for i = 1:n
    gatingGGIW{i}.r = zeros(0,1);
    gatingGGIW{i}.alpha = zeros(0,1);
    gatingGGIW{i}.beta = zeros(0,1);
    gatingGGIW{i}.x = zeros(x_dim,0);
    gatingGGIW{i}.P = zeros(x_dim,x_dim,0);
    gatingGGIW{i}.v = zeros(0,1);
    gatingGGIW{i}.V = zeros(d,d,0);
    for j = 1:plength
        if uniqueLabel(i) == label(j)
            gatingGroup{i} = union(gatingGroup{i},tentativeGroup{j});
            gatingGGIW{i}.r = [gatingGGIW{i}.r;ggiw_mb.r(j)];
            gatingGGIW{i}.alpha = [gatingGGIW{i}.alpha;ggiw_mb.alpha(j)];
            gatingGGIW{i}.beta = [gatingGGIW{i}.beta;ggiw_mb.beta(j)];
            gatingGGIW{i}.x = [gatingGGIW{i}.x ggiw_mb.x(:,j)];
            gatingGGIW{i}.P = cat(3,gatingGGIW{i}.P,ggiw_mb.P(:,:,j));
            gatingGGIW{i}.v = [gatingGGIW{i}.v;ggiw_mb.v(j)];
            gatingGGIW{i}.V = cat(3,gatingGGIW{i}.V,ggiw_mb.V(:,:,j));
        end
    end
    
end


end