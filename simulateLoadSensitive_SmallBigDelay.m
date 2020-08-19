function [average_Dep_Net_Congestion1,average_Dep_Net_Congestion2,average_Dep_Net_Congestion3...
    ,worst_Dep_Net_Congestion1,worst_Dep_Net_Congestion2,worst_Dep_Net_Congestion3...
    ] = simulateLoadSensitive_SmallBigDelay()

sim_num = 25;
duration = 3600;
deltaTime = 60;
iN = duration/deltaTime;

average_Dep_Net_Congestion1 = zeros(iN,1);
worst_Dep_Net_Congestion1 = zeros(iN,1);
average_Dep_Net_Uncertainty1 = zeros(iN+1,1);

average_Dep_Net_Congestion2 = zeros(iN,1);
worst_Dep_Net_Congestion2 = zeros(iN,1);
average_Dep_Net_Uncertainty2 = zeros(iN+1,1);

average_Dep_Net_Congestion3 = zeros(iN,1);
worst_Dep_Net_Congestion3 = zeros(iN,1);
% average_Dep_Net_Uncertainty3 = zeros(iN+1,1);

for i = 1:sim_num
    
    % Generate traffic
    N = 48;
    M = 140;
    D = duration;
    dT = deltaTime;
    tD = 10;
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
        AR = cat(3,5*ones(N,N,1),5*ones(N,N,D/dT-1));
    end
    
    DR = cat(3,zeros(N,N,floor(tD/dT)),dT./AR(:,:,1:D/dT-floor(tD/dT)));
    DAR = dT./AR-DR;
    DVAR = dT./AR+DR;
    DAR = cat(3,zeros(N,N,1),DAR,AR(:,:,D/dT-floor(tD/dT+1)));
    DVAR = cat(3,zeros(N,N,1),DVAR,AR(:,:,D/dT-floor(tD/dT+1)));
    DM = [3600*ones(5,iN);0.1*ones(N-5,iN)];
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
    % End of generate traffic
    
    
    %     Simulation type #1
    [Dep_Net_Congestion1,Dep_Net_Uncertainty1] = simulateLoadSensitive_1(iN,...
        capacities,ActualTraffic,PredictedTraffic,TrafficVariance3D,...
        TrafficVariance3D_Est,N,M);
    average_Dep_Net_Congestion1 = average_Dep_Net_Congestion1...
        + Dep_Net_Congestion1/sim_num;
    worst_Dep_Net_Congestion1 = max(worst_Dep_Net_Congestion1,Dep_Net_Congestion1);
    average_Dep_Net_Uncertainty1 = average_Dep_Net_Uncertainty1...
        + Dep_Net_Uncertainty1/sim_num;
    disp(['sim#1 ',num2str(i)]);
    
%     Simulation type #2
    [Dep_Net_Congestion2,Dep_Net_Uncertainty2] = simulateLoadSensitive_2(iN,...
        capacities,ActualTraffic,PredictedTraffic,TrafficVariance3D...
        ,TrafficVariance3D_Est,N,M);
    average_Dep_Net_Congestion2 = average_Dep_Net_Congestion2...
        + Dep_Net_Congestion2/sim_num;
    worst_Dep_Net_Congestion2 = max(worst_Dep_Net_Congestion2,Dep_Net_Congestion2);
    average_Dep_Net_Uncertainty2 = average_Dep_Net_Uncertainty2...
        + Dep_Net_Uncertainty2/sim_num;
    disp(['sim#2 ',num2str(i)]);
    
    File = fullfile(cd, ['Laptop Computer K',num2str(i)]);
    save(File);
    
    if i>1
        File_Old = fullfile(cd, ['Laptop Computer K',num2str(i-1),'.mat']);
        delete(File_Old);
    end

end

figure;plot(1:iN,average_Dep_Net_Congestion1);hold on;
plot(1:iN,average_Dep_Net_Congestion2);hold on;

figure;plot(1:iN,worst_Dep_Net_Congestion1);hold on;
plot(1:iN,worst_Dep_Net_Congestion2);hold on;

figure;plot(1:iN+1,average_Dep_Net_Uncertainty1);hold on;
plot(1:iN+1,average_Dep_Net_Uncertainty2);hold on;



end