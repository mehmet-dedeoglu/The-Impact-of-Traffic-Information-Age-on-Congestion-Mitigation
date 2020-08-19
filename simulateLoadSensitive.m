function res=simulateLoadSensitive()

% Inputs
% 
% Variables
% tD : Traffic duration
% DM : Delay mean matrix
% AR : Traffic arrival rate
% DAR : Difference of arrival rate
% Outputs
% 

N = 48;
M = 140;
D = 100;
dT = 2;
tD = 10;
iN = D/dT;
capacities = 1000*ones(M,1);

constant_simulate = 1;

% Generate Traffic
AR = cat(3,2*ones(N,N,1),2*ones(N,N,D/dT-1));
DR = cat(3,zeros(N,N,tD/dT),dT./AR(:,:,1:D/dT-tD/dT));
DAR = dT./AR-DR;
DVAR = dT./AR+DR;
DAR = cat(3,zeros(N,N,1),DAR,AR(:,:,D/dT-tD/dT+1));
DVAR = cat(3,zeros(N,N,1),DVAR,AR(:,:,D/dT-tD/dT+1));
DM = 10*ones(N,iN);
if constant_simulate == 0
    [TD,TO] = Inhomogenous_Poisson_Traffic(N,D,AR,dT,tD);
else
    [TD,TO] = Constant_Traffic(N,iN,dT);
    DAR = zeros(48,48,iN+2);
    DVAR = zeros(48,48,iN+2);
end
[RT,MST] = Measurement_Arrival_Times(N,1,DM,dT,iN);
[MeasuredTraffic,ActualTraffic,PredictedTraffic,TrafficVariance3D] ...
    = Traffic_Measurements(TD,TO,MST,iN,N,dT,DAR,DVAR);

% Generate Network
weightHist = [1*ones(3,M);zeros(iN-2,M)];
weightHist_Prop = weightHist;       % For proposed algorithm only
Can_Weights = weightHist_Prop(1,:); % For proposed algorithm only
choice = 1;

[~,G_ind,G_weight,~,G_adj,M] = createGraph(weightHist(1,:),choice,M);

% Find Shortest Paths
meanSPNumber = 2;
pathMatrix = zeros(meanSPNumber*N*N,N);
pathMatrixIndex = zeros(meanSPNumber*N*N,1);
pathCosts = zeros(meanSPNumber*N*N,1);
k = 1;
for j = 1:N
    for i = [1:j-1,j+1:N]
        [paths,pathCost] = dijkstra(G_ind,G_weight,j,i,N);
        [a,b] = size(paths);
        pathMatrix(k:k+a-1,:) = paths;
        pathMatrixIndex(k:k+a-1) = j;
        pathCosts(k:k+a-1) = pathCost;
        k = k + a;
    end
end

% Compute Routing Matrix
% PredictedTraffic(:,1) is only used to compute some metric we do not use
% here
[RM_Candidate,~,~,~,~,~,~] = ...
    flowDist(PredictedTraffic(:,1), pathMatrix, pathMatrixIndex,G_adj,pathCosts,N,M);
RM_Candidate_P = RM_Candidate; % For proposed algorithm only

% Reconfigure the network at the beginning
reconfig = zeros(iN+1,1);
reconfig(1) = 1;
reconfig_Iter = 1;
reconfig_P = zeros(iN+1,1);   % For proposed algorithm only
reconfig_P(1) = 1;            % For proposed algorithm only
reconfig_Iter_P = 1;          % For proposed algorithm only

% Simulation variables
Dep_Net_Cost = zeros(iN,1);
Dep_Net_Cost_P = zeros(iN,1); % For proposed algorithm only
flowSave = zeros(iN,M);         % For proposed algorithm only

% Controlling reconfigurations over iterations
for ell=1:iN
%     Update deployed routing matrix if reconfiguration decision is made in
%     the previous iteration #1
    if reconfig(ell)==1
        RM_Deployed = RM_Candidate;
    end
%     Compute corresponding costs related to deployed weights and candidate 
%     weights of previous iteration
%     Compute predicted link loads for the next iteration
    LL_Dep_Pred = RM_Deployed'*PredictedTraffic(:,ell);
    LL_Can_Pred = RM_Candidate'*PredictedTraffic(:,ell);
    
%     Compute link loads in the network for the current iteration
    LL_Dep_Net = RM_Deployed'*ActualTraffic(:,ell+1);
    LL_Can_Net = RM_Candidate'*ActualTraffic(:,ell+1);
    
%     Compute the latest measured link loads in the network
    LL_Dep_Meas = RM_Deployed'*MeasuredTraffic(:,ell+1);
    LL_Can_Meas = RM_Candidate'*MeasuredTraffic(:,ell+1);
    
%     Compute costs due to max link loads
    [pred_Dep_Cost,~] = max(LL_Dep_Pred.*LL_Dep_Pred...
        +diag(RM_Deployed'*TrafficVariance3D(:,:,ell)*RM_Deployed));
    [Dep_Net_Cost(ell),~] = max(LL_Dep_Net.*LL_Dep_Net);
    [pred_Can_Cost,~] = max(LL_Can_Pred.*LL_Can_Pred...
        +diag(RM_Candidate'*TrafficVariance3D(:,:,ell)*RM_Candidate));
    [Can_Net_Cost(ell),~] = max(LL_Can_Net.*LL_Can_Net); 
    
%     Compute a new set of candidate weights using LSAR (use link loads 
%     constituted if previous iteration's candidate weights are used)
    [weightHist]=LSAR(capacities,LL_Can_Pred,M,weightHist,ell+3);
    candidateWeights = weightHist(ell+3,:);
        
%     Compute corresponding costs related to candidate weights
%     Compute corresponding shortest paths
    [G,G_ind,G_weight,G_cap,G_adj,M] = createGraph(candidateWeights,choice,M);
    meanSPNumber = 2;
    pathMatrix = zeros(meanSPNumber*N*N,N);
    pathMatrixIndex = zeros(meanSPNumber*N*N,1);
    pathCosts = zeros(meanSPNumber*N*N,1);
    k = 1;
    for j = 1:N
        for i = [1:j-1,j+1:N]
            [paths,pathCost] = dijkstra(G_ind,G_weight,j,i,N);
            [a,b] = size(paths);
            pathMatrix(k:k+a-1,:) = paths;
            pathMatrixIndex(k:k+a-1) = j;
            pathCosts(k:k+a-1) = pathCost;
            k = k + a;
        end
    end
%     Compute corresponding routing table
    [RM_Candidate,~,~,~,~,~,~] = ...
    flowDist(PredictedTraffic(:,ell), pathMatrix, pathMatrixIndex,G_adj,pathCosts,N,M);

%     Compute predicted link loads for the next iteration
    LL_Can_Pred = RM_Candidate'*PredictedTraffic(:,ell);
    
%     Compute link loads in the network for the current iteration
    LL_Can_Net = RM_Candidate'*ActualTraffic(:,ell+1);
    
%     Compute the latest measured link loads in the network
    LL_Can_Meas = RM_Candidate'*MeasuredTraffic(:,ell+1);
    
%     Compute costs due to max link loads
    [pred_Can_Cost_No_Reconfig,~] = max(LL_Can_Pred.*LL_Can_Pred...
        +diag(RM_Candidate'*TrafficVariance3D(:,:,ell)*RM_Candidate));
    reconfig_Cost = 0.2/reconfig_Iter*pred_Dep_Cost;
    pred_Can_Cost = pred_Can_Cost_No_Reconfig + reconfig_Cost;
    reconfig_Iter = reconfig_Iter + 1;

%     Decide whether to reconfigure the network or not
    if pred_Dep_Cost > pred_Can_Cost
        reconfig(ell+1) = 1;
        reconfig_Iter = 1;
    end
    disp(ell)
    
    
    
    
    
%     Update deployed routing matrix if reconfiguration decision is made in
%     the previous iteration for proposed technique
    if reconfig_P(ell)==1
        RM_Deployed_P = RM_Candidate_P;
    end
%     Compute corresponding costs related to deployed weights
%     Compute predicted link loads for the next iteration
    LL_Dep_Pred_P = RM_Deployed_P'*PredictedTraffic(:,ell);
    LL_Can_Pred_P = RM_Candidate_P'*PredictedTraffic(:,ell);
    
%     Compute link loads in the network for the current iteration
    LL_Dep_Net_P = RM_Deployed_P'*ActualTraffic(:,ell+1);
    LL_Can_Net_P = RM_Candidate_P'*ActualTraffic(:,ell+1);
    
%     Compute the latest measured link loads in the network
    LL_Dep_Meas_P = RM_Deployed_P'*MeasuredTraffic(:,ell+1);
    LL_Can_Meas_P = RM_Candidate_P'*MeasuredTraffic(:,ell+1);
    
%     Compute costs due to max link loads
    [pred_Dep_Cost_P,~] = max(LL_Dep_Pred_P.*LL_Dep_Pred_P...
        +diag(RM_Deployed_P'*TrafficVariance3D(:,:,ell)*RM_Deployed_P));
    [Dep_Net_Cost_P(ell),~] = max(LL_Dep_Net_P.*LL_Dep_Net_P);
    [pred_Can_Cost_P,~] = max(LL_Can_Pred_P.*LL_Can_Pred_P...
        +diag(RM_Candidate_P'*TrafficVariance3D(:,:,ell)*RM_Candidate_P));
    [Can_Net_Cost_P(ell),~] = max(LL_Can_Net_P.*LL_Can_Net_P);
    
%     Compute a new set of candidate weights using proposed algorithm 
%  (use link loads constituted if previous iteration's candidate weights are used)
    [Can_Weights,flowSave]=...
        updateWeights_info(capacities,LL_Can_Pred_P,M,Can_Weights,flowSave);
    
%     Compute corresponding costs related to candidate weights
%     Compute corresponding shortest paths
    [G,G_ind,G_weight,G_cap,G_adj,M] = createGraph(Can_Weights,choice,M);
    meanSPNumber = 2;
    pathMatrix = zeros(meanSPNumber*N*N,N);
    pathMatrixIndex = zeros(meanSPNumber*N*N,1);
    pathCosts = zeros(meanSPNumber*N*N,1);
    k = 1;
    for j = 1:N
        for i = [1:j-1,j+1:N]
            [paths,pathCost] = dijkstra(G_ind,G_weight,j,i,N);
            [a,b] = size(paths);
            pathMatrix(k:k+a-1,:) = paths;
            pathMatrixIndex(k:k+a-1) = j;
            pathCosts(k:k+a-1) = pathCost;
            k = k + a;
        end
    end
%     Compute corresponding routing table
    [RM_Candidate_P,~,~,~,~,~,~] = ...
    flowDist(PredictedTraffic(:,ell), pathMatrix, pathMatrixIndex,G_adj,pathCosts,N,M);

%     Compute predicted link loads for the next iteration
    LL_Can_Pred_P = RM_Candidate_P'*PredictedTraffic(:,ell);
    
%     Compute link loads in the network for the current iteration
    LL_Can_Net_P = RM_Candidate_P'*ActualTraffic(:,ell+1);
    
%     Compute the latest measured link loads in the network
    LL_Can_Meas_P = RM_Candidate_P'*MeasuredTraffic(:,ell+1);
    
%     Compute costs due to max link loads
    [pred_Can_Cost_No_Reconfig_P,~] = max(LL_Can_Pred_P.*LL_Can_Pred_P...
        +diag(RM_Candidate_P'*TrafficVariance3D(:,:,ell)*RM_Candidate_P));
    reconfig_Cost_P = 0.2/reconfig_Iter_P*pred_Dep_Cost_P;
    pred_Can_Cost_P = pred_Can_Cost_No_Reconfig_P + reconfig_Cost_P;
    reconfig_Iter_P = reconfig_Iter_P + 1;

%     Decide whether to reconfigure the network or not
    if pred_Dep_Cost_P > pred_Can_Cost_P
        reconfig_P(ell+1) = 1;
        reconfig_Iter_P = 1;
    end
  
end

figure;plot(1:iN,Dep_Net_Cost);hold on;plot(1:iN,Dep_Net_Cost_P);







end