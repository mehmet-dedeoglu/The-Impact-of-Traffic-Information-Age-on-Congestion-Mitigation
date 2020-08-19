function [weights,flowSave]=...
    updateWeights_info(capacities,flows,M,weights,flowSave)

weightsOld = weights;

alpha = 4;
% weightSave = [weightSave(2:end,:);weights];
flowSave = [flowSave(2:end,:);flows'];
wmax = 100;
wmin = 4;
deltaWeights = zeros(1,M);
% std = sqrt(var(flows));
mean = sum(flows)/M;
mean2 = mean;
medW = median(weights);
cMean = capacities(1);
[fa,~] = size(flowSave);
fMeanVec = sum(flowSave(fa-4:end,:),1)/5;
for i = 1:100
    sPoint = ((alpha-1)*cMean-mean)/(alpha);
    flowTemp = flows(fMeanVec>sPoint);
    if sum(flowTemp) ~= 0
        alpha = 2*alpha;
        mean2 = sum(flowTemp)/length(flowTemp);
    end
end
deltaWeights(flows<mean) = -ceil(medW/5*(mean-flows(flows<mean))...
    /(mean));
deltaWeights(flows>sPoint) = ceil((flows(flows>sPoint)-...
    sPoint).*medW/(2*mean2));
% deltaWeights(flows>sPoint) = ceil((flows(flows>sPoint)-...
%     sPoint).*medW/(2*(cMean-sPoint)));
weights = weights + deltaWeights;
if weights == weightsOld
%     weights(flowsOld==max(flowsOld))=weights(flowsOld==max(flowsOld))+2;
    weights = weights - round(rand(1,M));
end

weights(weights>wmax) = wmax;
weights(weights<wmin) = wmin;

if max(weights)==wmax
    weights = ceil(weights/10);
end

end