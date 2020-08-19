function [Dep_Net_Congestion,Dep_Net_Uncertainty] = UnitOSPF(iN,...
        ActualTraffic,PredictedTraffic,TrafficVariance3D_Est,N,M)

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
[RM,~,~,~,~,~,~] = ...
    flowDist(PredictedTraffic(:,1), pathMatrix, pathMatrixIndex,G_adj,pathCosts,N,M);

% Simulation variables
Dep_Net_Congestion = zeros(iN,1);
Dep_Net_Uncertainty = zeros(iN+1,1);

% Controlling reconfigurations over iterations
for ell=1:iN
    
%     Compute link loads in the network for the current iteration
    LL_Dep_Net = RM'*ActualTraffic(:,ell+1);
    
%     Compute costs due to max link loads
    [Dep_Net_Congestion(ell),~] = max(LL_Dep_Net);

%     Compute actual variances for simulations
    [Dep_Net_Uncertainty(ell),~] = max(diag(RM'*...
        TrafficVariance3D_Est(:,:,ell+1)*RM));
    disp(ell)
  
end

% figure;plot(1:iN,Dep_Net_Cost);title('Deployed network maximum link load cost');drawnow;
% figure;plot(1:iN,Dep_Net_Congestion);title('Deployed network maximum link load');drawnow;

end