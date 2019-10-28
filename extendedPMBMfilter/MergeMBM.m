function [MBM_m,log_Wmbm_m] = MergeMBM(MBM,log_Wmbm,model)

% Code by Karl Granström

MBM_MergingThreshold = model.mergingThreshold;

% Number of MBs
Nmb = numel(MBM);

if Nmb<=1
    MBM_m = MBM;
    log_Wmbm_m = log_Wmbm;
else
    
    % Allocate memory
    kld = inf(Nmb,Nmb);
    thres = zeros(Nmb,Nmb);
    
    for imb = 1:Nmb-1
        % i:th MB
        MBi = MBM{imb};
        % Number of Bernoullis
        Nmbi = length(MBi.r);
        for jmb = imb+1:Nmb
            
            % j:th MB
            MBj = MBM{jmb};
            % Check if MBi and MBj have the same number of Bernoullis
            if Nmbi == length(MBj.r)
                % Compute the Kullback Leibler difference
                [kld(imb,jmb),~] = KLdiffMB(MBi,MBj);
                % Set the threshold proportional to the number of
                % Bernoullis
                thres(imb,jmb) = MBM_MergingThreshold*Nmbi;
            end
            kld(jmb,imb) = kld(imb,jmb);
            thres(jmb,imb) = thres(imb,jmb);
        end
    end
    
    % Neighbour matrix
    NeighbourMatrix = kld<thres;
    
    ptsC  = zeros(Nmb,1);
    C     = {};
    Nc    = 0;               % Cluster counter.
    Pvisit = zeros(Nmb,1);  % Array to keep track of points that have been visited.
    
    for n = 1:Nmb
        if ~Pvisit(n)                            % If this point not visited yet
            Pvisit(n) = 1;                       % mark as visited
            %neighbourPts = regionQuery(P, n, E); % and find its neighbours
            neighbourPts = find(NeighbourMatrix(n,:));
            
            Nc = Nc + 1;    % Increment number of clusters and process
            % neighbourhood.
            
            C{Nc} = [n];    % Initialise cluster Nc with point n
            ptsC(n) = Nc;   % and mark point n as being a member of cluster Nc.
            
            ind = 1;        % Initialise index into neighbourPts array.
            
            % For each point P' in neighbourPts ...
            while ind <= length(neighbourPts)
                
                nb = neighbourPts(ind);
                
                if ~Pvisit(nb)        % If this neighbour has not been visited
                    Pvisit(nb) = 1;   % mark it as visited.
                    
                    % Find the neighbours of this neighbour and if it has
                    % enough neighbours add them to the neighbourPts list
                    %neighbourPtsP = regionQuery(P, nb, E);
                    neighbourPtsP = find(NeighbourMatrix(nb,:));
                    neighbourPts = [neighbourPts  neighbourPtsP];
                end
                
                % If this neighbour nb not yet a member of any cluster add it
                % to this cluster.
                if ~ptsC(nb)
                    C{Nc} = [C{Nc} nb];
                    ptsC(nb) = Nc;
                end
                
                ind = ind + 1;  % Increment neighbour point index and process
                % next neighbour
            end
        end
    end
    
    % Allocate memory
    MBM_m = cell(numel(C),1);
    log_Wmbm_m = zeros(1,numel(C));
    
    % Iterate over clusters
    for ic = 1:numel(C)
        
        if length(C{ic})>1
            % If more than one in cluster
            
            % Find maximum weight MB
            [~,idxMax] = max(log_Wmbm(C{ic}));
            
            % Store this MB
            MBM_m{ic} = MBM{C{ic}(idxMax)};
            
            % Sum of weights
            [~,log_Wmbm_m(ic)] = normalizeLogWeights(log_Wmbm(C{ic}));
            
        else
            % Only one in cluster
            
            % Store this MB and weight
            MBM_m{ic} = MBM{C{ic}};
            log_Wmbm_m(ic) = log_Wmbm(C{ic});
        end
    end
    
    %log_Wmbm_m = normalizeLogWeights(log_Wmbm_m);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kld,assignment] = KLdiffMB(MBi,MBj)

KLdiff = zeros(length(MBi.r),length(MBj.r));

for ik = 1:length(MBi.r)
    for jk = 1:length(MBj.r)
        
        [KLdiff(ik,jk),~,~,~,~] = GGIW_KLdiff(...
            MBi.r(ik),MBi.alpha(ik),MBi.beta(ik),MBi.x(:,ik),MBi.P(:,:,ik),MBi.v(ik),MBi.V(:,:,ik),...
            MBj.r(jk),MBj.alpha(jk),MBj.beta(jk),MBj.x(:,jk),MBj.P(:,:,jk),MBj.v(jk),MBj.V(:,:,jk));
    end
end

assignment = auctionAlgorithm(-KLdiff);

kld = 0;

for ik = 1:length(assignment)
    kld = kld + KLdiff(ik,assignment(ik));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [KLdiff,w_diff,G_diff,N_diff,IW_diff] = ...
    GGIW_KLdiff(w1,a1,b1,m1,P1,v1,V1,w2,a2,b2,m2,P2,v2,V2)

[KLdiv_w1,KLdiv_G1,KLdiv_N1,KLdiv_IW1] = GGIW_KLdiv(w1,a1,b1,m1,P1,v1,V1,w2,a2,b2,m2,P2,v2,V2);
[KLdiv_w2,KLdiv_G2,KLdiv_N2,KLdiv_IW2] = GGIW_KLdiv(w2,a2,b2,m2,P2,v2,V2,w1,a1,b1,m1,P1,v1,V1);

w_diff = KLdiv_w1+KLdiv_w2;
G_diff = KLdiv_G1+KLdiv_G2;
N_diff = KLdiv_N1+KLdiv_N2;
IW_diff = KLdiv_IW1+KLdiv_IW2;

KLdiff = w_diff+G_diff+N_diff+IW_diff;
% KLdiff = N_diff+IW_diff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [KLdiv_w,KLdiv_G,KLdiv_N,KLdiv_IW] = GGIW_KLdiv(w1,a1,b1,m1,P1,v1,V1,w2,a2,b2,m2,P2,v2,V2)

% Function that computes the KL-divergence between two
% Gamma-Gaussian-inverse Wishart (GGIW) distributions
%
% p(g,x,X) = Gam(g;a1,b1)N(x;m1,P1)IW(X;v1,V1)
% q(g,x,X) = Gam(g;a2,b2)N(x;m2,P2)IW(X;v2,V2)
%
% D_KL = D(p||q) = int p ln(p/q) dx

% Dimension of extension
d = size(V1,1);
% Dimension of kinematical state
n_x = length(m1);

KLdiv_w = w1*(log(w1)-log(w2));

KLdiv_G = w1*(a1*log(b1)-a2*log(b2)+gammaln(a2)-gammaln(a1)...
    +(a1-a2)*(psi(0,a1)-log(b1))+a1*(b2/b1-1));

KLdiv_N = w1*(...
    -0.5*log(det(P1))+0.5*log(det(P2))...
    -0.5*n_x+0.5*(m1-m2)'*(P2\(m1-m2))...
    +0.5*trace(P2\P1)...
    );

KLdiv_IW = w1*(...
    0.5*(v1-d-1)*log(det(V1))-0.5*(v2-d-1)*log(det(V2))...
    +sum(gammaln((v2-d-(1:d))/2)-gammaln((v1-d-(1:d))/2))...
    +0.5*(v2-v1)*(log(det(V1))-sum(psi(0,(v1-d-(1:d))/2)))...
    +trace(-0.5*(v1-d-1)*(V1\(V1-V2)))...
    );