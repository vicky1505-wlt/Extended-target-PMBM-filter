function [ Error ]= OSPAmetric(X,Y,c,p)

% Compute OSPA distance between two finite sets X and Y
% X.x -- Nx * J matrix with kinematic vectors
% X.X -- d * d * J tensor with random matrices
%        c  -   cut-off parameter
%        p  -   p-parameter for the metric
% Output: scalar distance between X and Y
% Note: the Euclidean 2-norm is used as the "base" distance on the region


if isempty(X.x) && isempty(Y.x)
    dist = 0;
    
    Error = [dist,0,0];
    
    return;
end

if isempty(X.x) || isempty(Y.x)
    dist = c;
    
    Error = [dist,0,c];
    
    return;
end


%Calculate sizes of the input point patterns
n = size(X.x,2);
m = size(Y.x,2);

% Distance matrix
D = repmat(c,[n m]);

for ix = 1:n
    for iy = 1:m
        % Gaussian Wasserstein Distance for kinematics and extent
        gwd = GaussianWassersteinDistance(X.x(1:2,ix),X.X(:,:,ix),Y.x(1:2,iy),Y.X(:,:,iy));
        % Poisson rates
        %prd = abs(X.g(ix)-Y.g(iy));
        
        % Apply threshold c
        D(ix,iy) = gwd;
    end
end
D = min(c,D).^p;

%Compute optimal assignment and cost using the Hungarian algorithm
[~,cost]= Hungarian(D);

%Calculate final distance

location_error = (1/max(m,n)*cost)^(1/p);
cardinality_error = (1/max(m,n)*c^p*abs(m-n))^(1/p);

dist = location_error+cardinality_error;

Error = [dist,location_error,cardinality_error];

end
    