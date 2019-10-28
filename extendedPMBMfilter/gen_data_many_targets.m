function [groundTruth,Z,model] = gen_data_many_targets(para,targetTracks)

K = max([targetTracks.deathTime]);

% Generate model parameters, ground truth and measurements 
x_dim = 4;
z_dim = 2;
model.x_dim = x_dim;                           % dimension of target state
model.z_dim = z_dim;                           % dimension of measurements
d = 2;                                         % dimension of target extent
model.K = K;                                   % number of scans
model.Ps = 0.99;                               % target survival probability
model.Pd = para.detection_prob;

% motion model parameters for Gamma and inverse-Wishart distributions
model.eta = 2;
model.tao = 10;
% motion model parameters for Gaussian distribution
Ts = 1;
model.Ts = Ts;
model.F = [eye(2) Ts*eye(2);zeros(2) eye(2)];  % transition model
sigma_v = 0.1;                                 % process noise covariance         
model.Q = kron([Ts^4/4 Ts^3/2; Ts^3/2 Ts^2]*sigma_v^2,eye(2));

% measurement model
model.H = [1 0 0 0;0 1 0 0];
sigma_v = 0.1;
model.D= eye(z_dim)*sigma_v^2; 
model.R= model.D*model.D';                     % observation noise covariance

% target initial state birth/death time
nbirths = length(targetTracks);
xstart = zeros(x_dim,nbirths);
xstart(:,1) = [-75 -75 0 0];
xstart(:,2) = [-75 75 0 0];
xstart(:,3) = [75 75 0 0];
xstart(:,4) = [75 -75 0 0];

X = cell(K,1);
E = cell(K,1);
Gam = cell(K,1);
N = zeros(K,1);
groundTruth = cell(K,1);
% generate tracks (ground truth)
for targetnum = 1:nbirths
    for k = targetTracks(targetnum).birthTime:targetTracks(targetnum).deathTime
        targetstate = targetTracks(targetnum).x(1:4,k-targetTracks(targetnum).birthTime+1);
        targetextent = targetTracks(targetnum).X(:,:,k-targetTracks(targetnum).birthTime+1);
        X{k} = [X{k} targetstate];
        E{k} = cat(3,E{k},targetextent);
        Gam{k} = [Gam{k} targetTracks(targetnum).g];
        N(k) = N(k) + 1;
    end
end

for k = 1:K
    groundTruth{k}.x = X{k};
    groundTruth{k}.X = E{k};
end

range_c = [-1 1;-1 1]*200;              % surveillance area
lambda_c = para.false_alarm_rate;       % Poisson rate of clutter
model.lambda_fa= lambda_c/prod(range_c(:,2)-range_c(:,1));  % clutter intensity

% generate measurements
Z = cell(K,1);
for k = 1:K
    for i = 1:N(k)
        if rand < model.Pd              % simulate missed detection
            p = mvnrnd(X{k}(1:2,i),E{k}(:,:,i)+model.R,poissrnd(Gam{k}(i)));
            Z{k} = [Z{k} p'];
        end
    end
    % add false alarm
    N_c = poissrnd(lambda_c);
    C = repmat(range_c(:,1),[1 N_c])+ diag(range_c*[-1;1])*rand(2,N_c);
    Z{k} = [Z{k} C];
end

% initialise PPP birth intensity parameters
nbirths = 4;
model.alpha_b = 32*ones(nbirths,1);
model.beta_b = 4*ones(nbirths,1);
model.vb = 12*ones(nbirths,1);
model.wb = 0.01*ones(nbirths,1);

model.xb = zeros(x_dim,nbirths);
model.Pb = zeros(x_dim,x_dim,nbirths);
model.Vb = zeros(d,d,nbirths);
for i = 1:nbirths
    model.xb(:,i) = xstart(:,i);
    model.Pb(:,:,i) = diag([ 1; 1; 1; 1 ])*diag([ 1; 1; 1; 1 ])';
    model.Vb(:,:,i) = 12*eye(2);
end

% Thresholds
% Bernoullis with existence probability less than this threshold are pruned
model.threshold_r = 1e-3;
% PPP intensity components with weight less than this threshold are pruned
model.threshold_u = 1e-3;
model.threshold_recycle = 0.1;      % recycling threshold
% MBs that correspond to the likelihood 1-model.threshold_w are kept
model.threshold_w = 1e-3;
model.M = 100;                       % M-best data association, called by Murty
model.mergingThreshold = 0.05;       % merging threshold in MBM merging
% merging threshold for Bernoulli merging in the recycling step
model.mergingThreshold_recycle = 2; 

% Gating parameters
Pg = 0.999;                             %gate size in percentage
model.gamma= chi2inv(Pg,model.z_dim);   %gate size 

% db-scan parameters, a grid search for hyperparameters
model.max_dist = 5;
model.min_dist = 0.1;
model.grid_dist = 0.5;

end