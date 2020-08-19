function duration = MCF()

sim_num = 1;
duration = 3000;
deltaTime = 60;
iN = duration/deltaTime;

average_Dep_Net_Congestion1 = zeros(sim_num,1);
worst_Dep_Net_Congestion1 = zeros(sim_num,1);
average_Dep_Net_Uncertainty1 = zeros(sim_num,1);

average_Dep_Net_Congestion2 = zeros(sim_num,1);
worst_Dep_Net_Congestion2 = zeros(sim_num,1);
average_Dep_Net_Uncertainty2 = zeros(sim_num,1);

for ell = 1:sim_num

    criteria = 0;
    while criteria == 0
        [adjacency, ~] = waxtop(0.15,0.7,0.7,[0 10 0 10]);
        adjacency = full(adjacency);
        adjacency = round((adjacency + adjacency')/2);
        [N,~] = size(adjacency);%48
        Adj_temp = adjacency;
        Adj_temp(1:N+1:N*N) = 1;
        [~,~,rr,~] = dmperm(Adj_temp);
        if length(rr) == 2 && N == 13
            criteria = 1;
        end
    end

    % Generate traffic

        M = nnz(adjacency);%140
        D = duration;
        dT = deltaTime;
        tD = 1;
        iN = D/dT;
        capacities = 1000*ones(M,1);

        constant_simulate = 0;
        stationary = 1;

    %     Generate non-stationary traffic

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
        DM = 1000*ones(N,iN);
        if constant_simulate == 0
            [TD,TO] = Inhomogenous_Poisson_Traffic(N,D,AR,dT,tD);
        else
            [TD,TO] = Constant_Traffic(N,iN,dT);
            DAR = zeros(48,48,iN+2);
            DVAR = zeros(48,48,iN+2);
        end
        [~,MST] = Measurement_Arrival_Times(N,1,DM,dT,iN);
        [~,ActualTraffic,PredictedTraffic,TrafficVariance3D,TrafficVariance3D_Est] ...
            = Traffic_Measurements(TD,TO,MST,iN,N,dT,DAR,DVAR);
    %     End of generate traffic

    for i = iN-2

    %--------------------------------------------------------------------------
    %   Formulate optimization problem
    %     theta = zeros(M,1);
    %     gamma = zeros(N*N,1);
        A_B = find(adjacency');
        A_D = zeros(N,N);
        A_E = zeros(N,N);
        A_D_Len = zeros(N,1);
        A_E_Len = zeros(N,1);
    %     for s = 1:N
    %         for d = [1:s-1,s+1:N]
    %             gamma((s-1)*N+d) = 1;
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
    %                 theta(A_B==A_C) = 1;
    %                 theta(A_B==N*(0:N-1)+m) = -1;
                end
    %         end
    %     end
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
                    disp([s,d]);
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
                    disp([s,d]);
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
    %         (G2.*G2)*(PredictedTraffic(:,i+1).*PredictedTraffic(:,i+1)+...
    %             diag(TrafficVariance3D(:,:,i)))

            (G2*PredictedTraffic(:,i)).*(G2*PredictedTraffic(:,i))...
                +(G2.*G2)*diag(TrafficVariance3D(:,:,i))<=...
                capacities.*capacities*maxVal2;%
            G2 <= 1;
            G2 >= 0;
            maxVal2 >= 0;
            for bb = 1:N+1:N*N
                G2(:,bb) == 0;
            end
        cvx_end    

        Dep_Net_Congestion1 = max(G*ActualTraffic(:,i)./capacities);
        Dep_Net_Congestion2 = max(G2*ActualTraffic(:,i)./capacities);
        Dep_Net_Uncertainty1 = max((G.*G)*diag(TrafficVariance3D_Est(:,:,i)));
        Dep_Net_Uncertainty2 = max((G2.*G2)*diag(TrafficVariance3D_Est(:,:,i)));
        
            %     Simulation type #1
        average_Dep_Net_Congestion1 = average_Dep_Net_Congestion1...
            + Dep_Net_Congestion1/sim_num;
        worst_Dep_Net_Congestion1 = max(worst_Dep_Net_Congestion1,Dep_Net_Congestion1);
        average_Dep_Net_Uncertainty1 = average_Dep_Net_Uncertainty1...
            + Dep_Net_Uncertainty1/sim_num;
        disp(['sim#1 ',num2str(ell)]);

    %     Simulation type #2
        average_Dep_Net_Congestion2 = average_Dep_Net_Congestion2...
            + Dep_Net_Congestion2/sim_num;
        worst_Dep_Net_Congestion2 = max(worst_Dep_Net_Congestion2,Dep_Net_Congestion2);
        average_Dep_Net_Uncertainty2 = average_Dep_Net_Uncertainty2...
            + Dep_Net_Uncertainty2/sim_num;
        disp(['sim#2 ',num2str(ell)]);

%         File = fullfile(cd, ['Laptop Computer Tradeoff A',num2str(ell)]);
%         save(File);
% 
%         if ell>1
%             File_Old = fullfile(cd, ['Laptop Computer Tradeoff A',num2str(ell-1),'.mat']);
%             delete(File_Old);
%         end
        
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