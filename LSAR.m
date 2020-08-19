function [weightHist]=LSAR(capacities,linkFlows,M,weightHist,ell)

minWeightUpdate = 0.1; %0.2
maxWeight = 3;
b = 1;
a = 0.5;
weights = ones(1,M);
linkLoads = linkFlows./capacities;
weights(linkLoads>0.5) = (maxWeight-b)/(1-a)*(linkLoads(linkLoads>0.5)-a)+b;
% (c-b)/(1-a)*(x-a)+b=y

aver_num = 3;
weightHist(ell,:) = (weights+sum(weightHist(ell-aver_num:ell-1,:),1))/(aver_num+1);
bool = abs(weightHist(ell,:)-weightHist(ell-1,:))<weightHist(ell-1,:)*minWeightUpdate;
weightHist(ell,bool) = weightHist(ell-1,bool);

end