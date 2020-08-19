function [TD,TO]=Inhomogenous_Poisson_Traffic(N,D,AR,dT,tD)

% Inputs
% N         : (1x1)         : Number of routers in the network
% D         : (1x1)         : Total simulation duration
% AR        : (NxNxD/dT)    : Arrival rates
% dT        : (1x1)         : Duration of one iteration
% Variables
% i         : (1x1)         : Source index
% j         : (1x1)         : Destination index
% tChange   : (1x1)         : The time between current and next arrival
% k         : (1x1)         : Iteration that current arrival falls into
% n         : (1x1)         : Output index of current arrival
% maxRate   : (1x1)         : Maximum arrival rate of traffic demands
% t         : (1x1)         : One row of traffic arrival times
% Outputs
% TD        : (NxNx?)       : Traffic demand matrix

maxRate = max(max(max(AR)));
TD = zeros(N,N,round(2*D/maxRate));
TO = zeros(N,N,round(2*D/maxRate));
for i=1:N
    for j=1:N
        t = reshape(TD(i,j,:),[],1);
        k=1;
        n=1;
        flag=0;
        while flag<D
            n=n+1;
            tChange=exprnd(AR(i,j,k));
            if flag==0
                t(n) = t(n-1)+tChange;
            else
                t(n) = flag+tChange;
            end
            if t(n)>k*dT
                k=k+1;
                n=n-1;
                if flag==0
                    flag=dT*floor(t(n)/dT)+dT;
                else
                    flag=flag+dT;
                end
            else
                flag=0;
            end
        end
%         disp([i,j])
        [count,T_ordered]=Counting_Process(t,tD,D);
        TD(i,j,1:length(count))=count;
        TO(i,j,1:length(count))=T_ordered(1,:);
    end
end

end