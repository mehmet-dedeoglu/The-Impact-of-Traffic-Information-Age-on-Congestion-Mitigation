function [meanTrafficDemand,average_Dep_Net_Congestion1,average_Dep_Net_Congestion2] = LSSimulation()

sim_num = 1;
duration = 3600;
deltaTime = 30;
iN = duration/deltaTime;
delay = [0.01,100];
trafficNumber = length(delay);

average_Dep_Net_Congestion1 = zeros(iN,trafficNumber);
average_Dep_Net_Uncertainty1 = zeros(iN+1,trafficNumber);
reconfig1 = zeros(iN+1,trafficNumber);

average_Dep_Net_Congestion2 = zeros(iN,trafficNumber);
average_Dep_Net_Uncertainty2 = zeros(iN+1,trafficNumber);
reconfig2 = zeros(iN+1,trafficNumber);

average_Dep_Net_Congestion3 = zeros(iN,trafficNumber);
average_Dep_Net_Uncertainty3 = zeros(iN+1,trafficNumber);

average_Dep_Net_Congestion4 = zeros(iN,trafficNumber);
average_Dep_Net_Uncertainty4 = zeros(iN+1,trafficNumber);

average_Dep_Net_Congestion5 = zeros(iN,trafficNumber);
average_Dep_Net_Uncertainty5 = zeros(iN+1,trafficNumber);

average_Dep_Net_Congestion6 = zeros(iN,trafficNumber);
average_Dep_Net_Uncertainty6 = zeros(iN+1,trafficNumber);

average_Dep_Net_Congestion7 = zeros(iN,trafficNumber);
average_Dep_Net_Uncertainty7 = zeros(iN+1,trafficNumber);

for i = 1:sim_num
    
    % Generate traffic
    N = 48;%15;%14;%23;%48
    M = 140;%72;%46;%74;%140
    D = duration;
    dT = deltaTime;
    tD = 1;
    iN = D/dT;
    capacities = 1000*ones(M,1);

    constant_simulate = 1;
    stationary = 0;
    
%     Generate non-stationary traffic
    
    if stationary == 0
        unif_rand1 = 1*rand(N,N,1)+2;
        unif_rand2 = 1*rand(N,N,D/dT-1)+2;
        AR = cat(3,unif_rand1,unif_rand2);
    else
        AR = cat(3,0.01*ones(N,N,1),0.01*ones(N,N,D/dT-1));
    end
    
    DR = cat(3,zeros(N,N,floor(tD/dT)),...
        (ceil(tD/dT)*dT-tD)./AR(:,:,ones(1,ceil(tD/dT)-floor(tD/dT))),...
        dT./AR(:,:,1:D/dT-ceil(tD/dT)));
    DAR = dT./AR-DR;
    DVAR = dT./AR+DR;
    DAR = cat(3,zeros(N,N,1),DAR,AR(:,:,D/dT-floor(tD/dT+1)));
    DVAR = cat(3,zeros(N,N,1),DVAR,AR(:,:,D/dT-floor(tD/dT+1)));
    
    if constant_simulate == 0
        [TD,TO] = Inhomogenous_Poisson_Traffic(N,D,AR,dT,tD);
    else
        [TD,TO] = Constant_Traffic(N,iN,dT);
        DAR = zeros(N,N,iN+2);
        DVAR = zeros(N,N,iN+2);
    end
    
    for kk = 1:length(delay)
        
        DM = delay(kk)*ones(N,iN);
        [~,MST] = Measurement_Arrival_Times(N,1,DM,dT,iN);
        [~,ActualTraffic,PredictedTraffic,TrafficVariance3D,TrafficVariance3D_Est] ...
            = Traffic_Measurements(TD,TO,MST,iN,N,dT,DAR,DVAR);
        
        meanTrafficDemand = sum(ActualTraffic,1)/N/N;
        % End of generate traffic
    
        %     Simulation type #1 (simulate proposed system with average metric)
        [average_Dep_Net_Congestion1(:,kk),average_Dep_Net_Uncertainty1(:,kk)...
            ,reconfig1(:,kk)] = LoadSensitive_1(iN,capacities,ActualTraffic...
            ,PredictedTraffic,TrafficVariance3D,TrafficVariance3D_Est,N,M);
        disp(['sim#1 ',num2str(i)]);
    
%         %     Simulation type #2
%         [average_Dep_Net_Congestion2(:,kk),average_Dep_Net_Uncertainty2(:,kk)...
%             ,reconfig2(:,kk)] = LoadSensitive_2(iN,capacities,ActualTraffic...
%             ,PredictedTraffic,TrafficVariance3D,TrafficVariance3D_Est,N,M);
%         disp(['sim#2 ',num2str(i)]);
%         
        %     Simulation UnitOSPF
        [average_Dep_Net_Congestion3(:,kk),average_Dep_Net_Uncertainty3(:,kk)] = UnitOSPF(iN,...
            ActualTraffic,PredictedTraffic,TrafficVariance3D_Est,N,M);
        disp(['sim UnitOSPF ',num2str(i)]);
%         
%         %     Simulation MCF full mean
%         [average_Dep_Net_Congestion4(:,kk),average_Dep_Net_Uncertainty4(:,kk)] = MCFFullMean(iN,...
%             capacities,ActualTraffic,PredictedTraffic,TrafficVariance3D...
%             ,TrafficVariance3D_Est,N,M);
%         disp(['sim# MCF Full Mean',num2str(i)]);
%         
        %     Simulation MCF partial mean
        [average_Dep_Net_Congestion5(:,kk),average_Dep_Net_Uncertainty5(:,kk)] = MCFPartialMean(iN,...
            capacities,ActualTraffic,PredictedTraffic,TrafficVariance3D...
            ,TrafficVariance3D_Est,N,M);
        disp(['sim# MCF Partial Mean',num2str(i)]);
%         
%         %     Simulation MCF full variance
%         [average_Dep_Net_Congestion6(:,kk),average_Dep_Net_Uncertainty6(:,kk)] = MCFFullVar(iN,...
%             capacities,ActualTraffic,PredictedTraffic,TrafficVariance3D...
%             ,TrafficVariance3D_Est,N,M);
%         disp(['sim# MCF Full Mean',num2str(i)]);
%         
%         %     Simulation MCF partial variance
%         [average_Dep_Net_Congestion7(:,kk),average_Dep_Net_Uncertainty7(:,kk)] = MCFPartialVar(iN,...
%             capacities,ActualTraffic,PredictedTraffic,TrafficVariance3D...
%             ,TrafficVariance3D_Est,N,M);
%         disp(['sim# MCF Partial Mean',num2str(i)]);
%         
%         
%         File = fullfile(cd, ['Desktop Simulation 2 A',num2str(i)]);
%         save(File);
%         
%         if i>1
%             File_Old = fullfile(cd, ['Desktop Simulation 2 A',num2str(i-1),'.mat']);
%             delete(File_Old);
%         end
    end
end

update1 = [average_Dep_Net_Congestion1(1,1);average_Dep_Net_Congestion1(reconfig1(2:end,1)==1,1)];
update2 = [average_Dep_Net_Congestion1(1,2);average_Dep_Net_Congestion1(reconfig1(2:end,2)==1,2)];
tt = 1:iN;
tt1 = [tt(1),tt(reconfig1(2:end,1)==1)];
tt2 = [tt(1),tt(reconfig1(2:end,2)==1)];

figure;plot(tt,average_Dep_Net_Congestion1(:,1),'LineWidth',1.5);hold on;
plot(1:iN,average_Dep_Net_Congestion1(:,2),'LineWidth',1.5);hold on;
plot(1:iN,average_Dep_Net_Congestion3(:,1),'LineWidth',1.5);hold on;
plot(tt1, update1,'b*','MarkerSize',10,'LineWidth',1.5);hold on;
plot(tt2, update2,'r+','MarkerSize',10,'LineWidth',1.5);grid on;
xlabel('Iteration');ylabel('Maximum Link Load');

% figure;
% plot(1:iN,reconfig1(1:iN,1));hold on;
% plot(1:iN,reconfig1(1:iN,2));hold on;



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

% name = 'Office Computer A';
% save(name);

end