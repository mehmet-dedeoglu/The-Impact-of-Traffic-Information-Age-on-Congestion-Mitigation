function [RT,MST] = Measurement_Arrival_Times(N,PI,DM,dT,IterNum)

% Input Parameters:
% N     : (1x1) : Number of nodes in the network.
% PI    : (Nx1) : Policy index vector that determines the policy for each
%                 node.
% DM    : (NxIterNum) : Delay matrix consisting of mean delay values.
% 
% Output Parameters:
% RT    : (Nx1) : Receive times of measurements at the controller.
if PI == 1
%   WT  : (Nx1) : Waiting time policy for each node.
    WT = zeros(N,1);
    
end

% Add communication delay to the controller.

RT = zeros(N,IterNum);
MST = zeros(N,IterNum);
for i = 1:IterNum
    UTT = (i-1)*dT*ones(N,1);
    UTTOld = UTT;
    tUp = i*dT;
    UTTN = zeros(N,1);
    index = true(N,1);
    while min(UTTN) < tUp
        RD = exprnd(DM(:,i));
        UTTN(index) = UTT(index)+ WT(index) + RD(index);
        index = UTTN < tUp;
        UTTOld(index) = UTT(index);
        UTT(index) = UTTN(index);
    end
    RT(:,i) = UTT;
    MST(:,i) = UTTOld;
    if i>1
        MST(RT(:,i)==(i-1)*dT,i) = MST(RT(:,i)==(i-1)*dT,i-1);
        MST(MST(:,i)==(i-1)*dT,i) = RT(MST(:,i)==(i-1)*dT,i-1);
        RT(RT(:,i)==(i-1)*dT,i) = RT(RT(:,i)==(i-1)*dT,i-1);
    end
end
RT = [zeros(N,1),RT];
MST = [zeros(N,1),MST];

end