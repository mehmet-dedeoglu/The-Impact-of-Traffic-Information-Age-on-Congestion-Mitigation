function [Dep_Net_Congestion,Dep_Net_Uncertainty] = MCFPartialVar(iN,...
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
RM_Deployed = RM_Candidate;

% Simulation variables
Dep_Net_Congestion = zeros(iN,1);
Dep_Net_Uncertainty = zeros(iN+1,1);

reconfigTimes = [2,12];

% Controlling reconfigurations over iterations
for ell=1:iN
%     Update deployed routing matrix if reconfiguration decision is made in
%     the previous iteration #1
    if ell == reconfigTimes(1) || ell == reconfigTimes(2)
        RM_Deployed = RM_Candidate;
    end
    
%     Compute link loads in the network for the current iteration
    LL_Dep_Net = RM_Deployed'*ActualTraffic(:,ell+1);
    
%     Compute costs due to max link loads
    [Dep_Net_Congestion(ell),~] = max(LL_Dep_Net);

%     Compute actual variances for simulations
    [Dep_Net_Uncertainty(ell),~] = max(diag(RM_Deployed'*...
        TrafficVariance3D_Est(:,:,ell+1)*RM_Deployed));

%     Compute a new set of candidate weights using LSAR (use link loads 
%     constituted if previous iteration's candidate weights are used)
%--------------------------------------------------------------------------
    %   Formulate optimization problem
    if ell == (reconfigTimes(1)-1) || ell == (reconfigTimes(2)-1)
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
        %     Solve optimization problem using CVX for second order moment metric
        diagVar = diag(TrafficVariance3D(:,:,ell));
        cvx_begin
        cvx_precision low
        cvx_solver SDPT3
        variable G2(M,N*N)
        variable maxVal2
        minimize maxVal2
        subject to
        for s = 1:N
            for d = [1:s-1,s+1:N]
                for cc = 1:N
                    if cc==s
                        sum(G2(A_D(cc,1:A_D_Len(cc)),(s-1)*N+d))-sum(G2(...
                            A_E(cc,1:A_E_Len(cc)),(s-1)*N+d))==1;
                    elseif cc~=d
                        sum(G2(A_D(cc,1:A_D_Len(cc)),(s-1)*N+d))-sum(G2(...
                            A_E(cc,1:A_E_Len(cc)),(s-1)*N+d))==0;
                    end
                end
            end
        end
        (G2*PredictedTraffic(:,ell)).*(G2*PredictedTraffic(:,ell))...
            +(G2.*G2)*diagVar<=capacities.*capacities*maxVal2;%
        G2 <= 1;
        G2 >= 0;
        maxVal2 >= 0;
        for bb = 1:N+1:N*N
            G2(:,bb) == 0;
        end
        cvx_end
        
        %     Compute corresponding routing table
        RM_Candidate = G2';
        vars1 = {'G2'};
        clear(vars1{:})
    end
    
    disp(ell)
  
end

% figure;plot(1:iN,Dep_Net_Cost);title('Deployed network maximum link load cost');drawnow;
% figure;plot(1:iN,Dep_Net_Congestion);title('Deployed network maximum link load');drawnow;

end