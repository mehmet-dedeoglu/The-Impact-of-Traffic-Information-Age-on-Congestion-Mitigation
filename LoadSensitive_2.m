function [Dep_Net_Congestion,Dep_Net_Uncertainty,reconfig] = LoadSensitive_2(iN,...
        capacities,ActualTraffic,PredictedTraffic,TrafficVariance3D,TrafficVariance3D_Est,N,M)

% This function computes reconfiguration decisions based on link load plus
% uncertainty costs. Further, new link weights are computed using the
% combination of link loads and uncertainty costs as well.

% Inputs
% 
% Variables
% tD : Traffic duration
% DM : Delay mean matrix
% AR : Traffic arrival rate
% DAR : Difference of arrival rate
% Outputs
% 

% Generate Network
weightHist = [1*ones(3,M);zeros(iN-2,M)];
choice = 4;

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
        [a,~] = size(paths);
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

% Reconfigure the network at the beginning
reconfig = zeros(iN+1,1);
reconfig(1) = 1;
reconfig_Iter = 1;

% Simulation variables
Dep_Net_Cost = zeros(iN,1);
Dep_Net_Congestion = zeros(iN,1);
Dep_Net_Uncertainty = zeros(iN+1,1);

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
%     LL_Can_Net = RM_Candidate'*ActualTraffic(:,ell+1);
    
% %     Compute the latest measured link loads in the network
%     LL_Dep_Meas = RM_Deployed'*MeasuredTraffic(:,ell+1);
%     LL_Can_Meas = RM_Candidate'*MeasuredTraffic(:,ell+1);
    
%     Compute costs due to max link loads
    [pred_Dep_Cost,~] = max(LL_Dep_Pred.*LL_Dep_Pred...
        +diag(RM_Deployed'*TrafficVariance3D(:,:,ell)*RM_Deployed));
    [Dep_Net_Cost(ell),~] = max(LL_Dep_Net.*LL_Dep_Net);
    [Dep_Net_Congestion(ell),~] = max(LL_Dep_Net);
%     [pred_Can_Cost,~] = max(LL_Can_Pred.*LL_Can_Pred...
%         +diag(RM_Candidate'*TrafficVariance3D(:,:,ell)*RM_Candidate));
%     [Can_Net_Cost(ell),~] = max(LL_Can_Net.*LL_Can_Net); 

%     Compute actual variances for simulations
    [Dep_Net_Uncertainty(ell),~] = max(diag(RM_Deployed'*...
        TrafficVariance3D_Est(:,:,ell+1)*RM_Deployed));
    
%     First compute link load+uncertaint cost
    track_load = LL_Can_Pred.*LL_Can_Pred;
    track_uncertainty = diag(RM_Candidate'*TrafficVariance3D(:,:,ell)*RM_Candidate);
    pred_Can_WUC = sqrt(LL_Can_Pred.*LL_Can_Pred...
        +diag(RM_Candidate'*TrafficVariance3D(:,:,ell)*RM_Candidate));

%     Compute a new set of candidate weights using LSAR (use link loads 
%     constituted if previous iteration's candidate weights are used)
    [weightHist]=LSAR(capacities,pred_Can_WUC,M,weightHist,ell+3);
    candidateWeights = weightHist(ell+3,:);
        
%     Compute corresponding costs related to candidate weights
%     Compute corresponding shortest paths
    [~,G_ind,G_weight,~,G_adj,M] = createGraph(candidateWeights,choice,M);
    meanSPNumber = 2;
    pathMatrix = zeros(meanSPNumber*N*N,N);
    pathMatrixIndex = zeros(meanSPNumber*N*N,1);
    pathCosts = zeros(meanSPNumber*N*N,1);
    k = 1;
    for j = 1:N
        for i = [1:j-1,j+1:N]
            [paths,pathCost] = dijkstra(G_ind,G_weight,j,i,N);
            [a,~] = size(paths);
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
    
% %     Compute link loads in the network for the current iteration
%     LL_Can_Net = RM_Candidate'*ActualTraffic(:,ell+1);
    
% %     Compute the latest measured link loads in the network
%     LL_Can_Meas = RM_Candidate'*MeasuredTraffic(:,ell+1);
    
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
  
end

% figure;plot(1:iN,Dep_Net_Cost);title('Deployed network maximum link load cost');drawnow;
% figure;plot(1:iN,Dep_Net_Congestion);title('Deployed network maximum link load');drawnow;

end