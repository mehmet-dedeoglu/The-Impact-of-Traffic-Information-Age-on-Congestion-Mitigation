function [Dep_Net_Congestion,Dep_Net_Uncertainty] = MCFFullMean(iN,...
        capacities,ActualTraffic,PredictedTraffic,TrafficVariance3D,TrafficVariance3D_Est,N,M)

% This function computes reconfiguration decisions based on just link loads. 
% Further, new link weights are computed using the link loads as well.

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
adjacency = G_adj;

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
    
%     Compute link loads in the network for the current iteration
    LL_Dep_Net = RM_Deployed'*ActualTraffic(:,ell+1);
    
%     Compute costs due to max link loads
    [pred_Dep_Cost,~] = max(LL_Dep_Pred.*LL_Dep_Pred);
    [Dep_Net_Cost(ell),~] = max(LL_Dep_Net.*LL_Dep_Net);
    [Dep_Net_Congestion(ell),~] = max(LL_Dep_Net);

%     Compute actual variances for simulations
    [Dep_Net_Uncertainty(ell),~] = max(diag(RM_Deployed'*...
        TrafficVariance3D_Est(:,:,ell+1)*RM_Deployed));

%     Compute a new set of candidate weights using LSAR (use link loads 
%     constituted if previous iteration's candidate weights are used)
%--------------------------------------------------------------------------
    %   Formulate optimization problem
    A_B = find(adjacency');
    A_D = zeros(N,N);
    A_E = zeros(N,N);
    A_D_Len = zeros(N,1);
    A_E_Len = zeros(N,1);
    for m = 1:N
        aa = 1;
        A_C = A_B(N*(m-1)+1<=A_B&A_B<=N*m);
        A_D_Len(m) = length(A_C);
        for mm = 1:A_D_Len(m)
            A_D(m,mm) = find(A_B==A_C(mm));
        end
        for mmm = 1:N
            temp = find(A_B==N*(mmm-1)+m);
            if isempty(temp)==0
                A_E(m,aa) = temp;
                aa = aa + 1;
            end
        end
        A_E_Len(m) = aa-1;
    end
    %--------------------------------------------------------------------------
    %     Solve optimization problem using CVX for first order moment metric
    cvx_begin
    cvx_precision low
    cvx_solver SDPT3
    variable G(M,N*N)
    variable maxVal
    minimize maxVal
    subject to
    for s = 1:N
        disp(s)
        for d = [1:s-1,s+1:N]
            %                         disp([s,d]);
            for cc = 1:N
                if cc==s
                    sum(G(A_D(cc,1:A_D_Len(cc)),(s-1)*N+d))-sum(G(...
                        A_E(cc,1:A_E_Len(cc)),(s-1)*N+d))==1;
                elseif cc~=d
                    sum(G(A_D(cc,1:A_D_Len(cc)),(s-1)*N+d))-sum(G(...
                        A_E(cc,1:A_E_Len(cc)),(s-1)*N+d))==0;
                end
            end
        end
    end
    G*PredictedTraffic(:,ell) <= capacities*maxVal;
    G <= 1;
    G >= 0;
    maxVal >= 0;
    for bb = 1:N+1:N*N
        G(:,bb) == 0;
    end
    cvx_end

%     Compute corresponding routing table
    RM_Candidate = G';
    vars1 = {'G'};
    clear(vars1{:})

%     Compute predicted link loads for the next iteration
    LL_Can_Pred = RM_Candidate'*PredictedTraffic(:,ell);
    
%     Compute costs due to max link loads
    [pred_Can_Cost_No_Reconfig,~] = max(LL_Can_Pred.*LL_Can_Pred);
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