function [MeasuredTraffic,ActualTraffic,PredictedTraffic,TrafficVariance3D,TrafficVariance3D_Est]...
    = Traffic_Measurements(TD,TO,MST,iN,N,dT,DAR,DVAR)

% Inputs
% TD        : (NxNx?)       : Traffic demand matrix
% MST       : ()            : Measurement sample times
% Variables
% 
% Outputs
% 

age_constant = ones(N,1)*(0:iN);
age_constant = dT*age_constant;
age = age_constant-MST;

MeasuredTraffic = zeros(N*N,iN+1);
ActualTraffic = zeros(N*N,iN+1);
PredictedTraffic = zeros(N*N,iN+1);

TrafficVariance = zeros(N*N,iN+1);
TrafficVarianceEst = zeros(N*N,iN+1);
TrafficVariance3D = zeros(N*N,N*N,iN+1);
TrafficVariance3D_Est = zeros(N*N,N*N,iN+1);

for i = 1:N
    disp(i)
    for j = 1:N
        TO_vec = reshape(TO(i,j,:),[],1);
        TD_vec = reshape(TD(i,j,:),[],1);
        for ell = 1:iN+1
            [~,ind1] = max(TO_vec(TO_vec<MST(i,ell)));
            MeasuredTraffic((i-1)*N+j,ell) = max([0,TD_vec(ind1)]);
            PredictedTraffic((i-1)*N+j,ell) = MeasuredTraffic((i-1)*N+j,ell)...
                + Mean_Change(age(i,ell),dT,reshape(DAR(i,j,:),[],1),ell);
            [~,ind2] = max(TO_vec(TO_vec<age_constant(i,ell)));
            ActualTraffic((i-1)*N+j,ell) = max([0,TD_vec(ind2)]);
            
            [TrafficVariance((i-1)*N+j,ell),TrafficVarianceEst((i-1)*N+j,ell)] = ...
                Variance_Change(age(i,ell),dT,reshape(DVAR(i,j,:),[],1),ell);
        end
    end
end
for aa=1:iN+1
    TrafficVariance3D(:,:,aa) = diag(TrafficVariance(:,aa));
    TrafficVariance3D_Est(:,:,aa) = diag(TrafficVarianceEst(:,aa));
end

end