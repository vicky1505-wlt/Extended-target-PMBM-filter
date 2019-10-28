function [ est ] = state_extract_MAP_MBM( W, ggiw_mbm )

if isempty(W)
    est.x = zeros(4,0);
    est.X = zeros(2,2,0);
    return;
end

% d = 2;
% extract number of MB components
num = length(W);
pcard = cell(num,1);
len = zeros(num,1);
max_card = zeros(num,1);

% for each MB component, obtain cardinality pmf
for i = 1:num
    r = ggiw_mbm{i}.r;
    r(r>1-eps) = 1-eps;
    pcard{i} = prod(1-r)*poly(-r./(1-r));
    len(i) = length(pcard{i});
    [~,n] = max(pcard{i});
    max_card(i) = n-1;
end

% add 0 to make each MB component has equal length of cardinality pmf
len_max = max(len);
pcard_mb = zeros(num,len_max);
for i = 1:num
    if len(i) < len_max
        pcard_mb(i,:) = W(i)*[pcard{i} zeros(1,len_max-length(pcard{i}))];
    else
        pcard_mb(i,:) = W(i)*pcard{i};
    end
end

% the MBM cardinality pmf is equal to the sum of the MB weights multiplied 
% with the MB cardinality PMFs
pcard_mbm = sum(pcard_mb);

% MAP cardinality estimate
[~,n] = max(pcard_mbm);
C_max = n-1;

% for each MB hypothesis in the MBM, find the MAP cardinality estimate
C_mb = zeros(num,1);
for i = 1:num
%     C_mb(i) = sum(ggiw_mbm{i}.r>=0.5);
    C_mb(i) = max_card(i);
end

% take the highest weight MB component satisfy C_max = C_mb
W(C_mb~=C_max) = 0;
[~,idx] = max(W);

% state extraction
est = state_extract_MAP(ggiw_mbm{idx});

end