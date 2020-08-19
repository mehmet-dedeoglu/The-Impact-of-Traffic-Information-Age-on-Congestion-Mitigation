function duration = MCFtemp()

nodeNumber = 10;
graphNumber = 25;
sampleNumber = 100;
duration = 1830;
deltaTime = 30;
delay = [0.01,0.1,1,10,100,1000];
% delay = [10,1000];
trafficNumber = length(delay);

average_Dep_Net_Congestion1 = zeros(graphNumber,length(delay));
worst_Dep_Net_Congestion1 = zeros(graphNumber,length(delay));
average_Dep_Net_Uncertainty1 = zeros(graphNumber,length(delay));

average_Dep_Net_Congestion2 = zeros(graphNumber,length(delay));
worst_Dep_Net_Congestion2 = zeros(graphNumber,length(delay));
average_Dep_Net_Uncertainty2 = zeros(graphNumber,length(delay));

for ii = 1:graphNumber
%     Generate a random graph
    criteria = 0;
    while criteria == 0
        [adjacency, ~] = waxtop(0.12,0.6,0.6,[0 10 0 10]);
        adjacency = full(adjacency);
        adjacency = round((adjacency + adjacency')/2);
        [N,~] = size(adjacency);%48
        Adj_temp = adjacency;
        Adj_temp(1:N+1:N*N) = 1;
        [~,~,rr,~] = dmperm(Adj_temp);
        if length(rr) == 2 && N == nodeNumber
            criteria = 1;
        end
    end
    
    for jj = 1:sampleNumber
        
        M = nnz(adjacency);%140
        D = duration;
        dT = deltaTime;
        tD = 1;
        iN = D/dT;
        capacities = 1000*ones(M,1);
        
        constant_simulate = 0;
        stationary = 1;
        
        % Generate non-stationary traffic
        
        if stationary == 0
            unif_rand1 = 9.9*rand(N,N,1)+0.1;
            unif_rand2 = 9.9*rand(N,N,D/dT-1)+0.1;
            AR = cat(3,unif_rand1,unif_rand2);
        else
            AR = cat(3,0.01*ones(N,N,1),0.01*ones(N,N,D/dT-1));
        end
        
        DR = cat(3,zeros(N,N,floor(tD/dT)),...
            (ceil(tD/dT)*dT-tD)./AR(:,:,ones(1,ceil(tD/dT)-floor(tD/dT))),...
            dT./AR(:,:,1:D/dT-ceil(tD/dT)));
        DAR = dT./AR-DR;
        DVAR = dT./AR+DR;
        DAR = cat(3,zeros(N,N,1),DAR,AR(:,:,D/dT-floor(tD/dT+1)));%????
        DVAR = cat(3,zeros(N,N,1),DVAR,AR(:,:,D/dT-floor(tD/dT+1)));%????
        
        % Generate traffic sample
        if constant_simulate == 0
            [TD,TO] = Inhomogenous_Poisson_Traffic(N,D,AR,dT,tD);
        else
            [TD,TO] = Constant_Traffic(N,iN,dT);
            DAR = zeros(48,48,iN+2);
            DVAR = zeros(48,48,iN+2);
        end
        
        for kk = 1:trafficNumber

            DM = delay(kk)*ones(N,iN);
            
            [~,MST] = Measurement_Arrival_Times(N,1,DM,dT,iN);
            [~,ActualTraffic,PredictedTraffic,TrafficVariance3D,TrafficVariance3D_Est] ...
                = Traffic_Measurements(TD,TO,MST,iN,N,dT,DAR,DVAR);

            %     End of generate traffic

            for i = iN-2
                diagVar = diag(TrafficVariance3D(:,:,i));
                diagVarEst = diag(TrafficVariance3D_Est(:,:,i));

    %             Delete useless variables
                vars = {'TrafficVariance3D','TrafficVariance3D_Est'};
                clear(vars{:})

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
                    G*PredictedTraffic(:,i) <= capacities*maxVal;
                    G <= 1;
                    G >= 0;
                    maxVal >= 0;
                    for bb = 1:N+1:N*N
                        G(:,bb) == 0;
                    end
                cvx_end    

                Dep_Net_Congestion1 = max(G*ActualTraffic(:,i));
                Dep_Net_Uncertainty1 = max((G.*G)*diagVarEst);
                vars1 = {'G'};
                clear(vars1{:})

            % -------------------------------------------------------------------------
            %     Solve optimization problem using CVX for second order moment metric
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
                    (G2*PredictedTraffic(:,i)).*(G2*PredictedTraffic(:,i))...
                        +(G2.*G2)*diagVar<=capacities.*capacities*maxVal2;%
                    G2 <= 1;
                    G2 >= 0;
                    maxVal2 >= 0;
                    for bb = 1:N+1:N*N
                        G2(:,bb) == 0;
                    end
                cvx_end    

                Dep_Net_Congestion2 = max(G2*ActualTraffic(:,i));
                Dep_Net_Uncertainty2 = max((G2.*G2)*diagVarEst);
                vars2 = {'G2'};
                clear(vars2{:})
                
            %     Simulation type #1
                average_Dep_Net_Congestion1(ii,kk) = average_Dep_Net_Congestion1(ii,kk)...
                    + Dep_Net_Congestion1/sampleNumber;
                worst_Dep_Net_Congestion1(ii,kk) = max(worst_Dep_Net_Congestion1(ii,kk),Dep_Net_Congestion1);
                average_Dep_Net_Uncertainty1(ii,kk) = average_Dep_Net_Uncertainty1(ii,kk)...
                    + Dep_Net_Uncertainty1/sampleNumber;
                disp(['sim# ',num2str(ii),'--',num2str(kk),'--',num2str(jj)]);

            %     Simulation type #2
                average_Dep_Net_Congestion2(ii,kk) = average_Dep_Net_Congestion2(ii,kk)...
                    + Dep_Net_Congestion2/sampleNumber;
                worst_Dep_Net_Congestion2(ii,kk) = max(worst_Dep_Net_Congestion2(ii,kk),Dep_Net_Congestion2);
                average_Dep_Net_Uncertainty2(ii,kk) = average_Dep_Net_Uncertainty2(ii,kk)...
                    + Dep_Net_Uncertainty2/sampleNumber;

                File = fullfile(cd, ['Desktop Computer Tradeoff A',num2str(ii),...
                    '--',num2str(kk),'--',num2str(jj)]);
                save(File);

                if jj>1
                    File_Old = fullfile(cd, ['Desktop Computer Tradeoff A',num2str(ii),...
                    '--',num2str(kk),'--',num2str(jj-1),'.mat']);
                    delete(File_Old);
                end
            end
        end
    end
end

% figure;plot(1:iN,average_Dep_Net_Congestion1);hold on;
% plot(1:iN,average_Dep_Net_Congestion2);hold on;
% % plot(1:iN,average_Dep_Net_Congestion3);
% 
% figure;plot(1:iN,worst_Dep_Net_Congestion1);hold on;
% plot(1:iN,worst_Dep_Net_Congestion2);hold on;
% % plot(1:iN,worst_Dep_Net_Congestion3);
% 
% figure;plot(1:iN+1,average_Dep_Net_Uncertainty1);hold on;
% plot(1:iN+1,average_Dep_Net_Uncertainty2);hold on;
% % plot(1:iN,average_Dep_Net_Uncertainty3);

end