function [routingMatrix,meanNumberOfSP,overallMeanNumberOfSP,...
    meanNumberOfSPW,overallMeanNumberOfSPW,meanSPlinkNumber,overallMeanSPlinkNumber] = ...
    flowDist(demand, pathMatrix, pathMatrixIndex,G_adj,pathCosts,N,M)

links = find(G_adj);
routingMatrix = zeros(N*N,M);
% % flows = zeros(N,N);
% numberOfSP = zeros(M,N*round(max(max(pathCosts))));
% numberOfSPW = zeros(M,N*round(max(max(pathCosts))));
% SPlinkNumber = zeros(M,N*round(max(max(pathCosts))));
for j = 1:N
    currentSource = pathMatrix(pathMatrixIndex == j,:);
%     costSource = pathCosts(pathMatrixIndex == j,:);
    for i = [1:j-1,j+1:N]
        currentDest = currentSource(currentSource(:,1) == i,:);
%         cost = costSource(currentSource(:,1) == i);
        [a,~] = size(currentDest);
        for k = 1:a
%             flowsUpdate = zeros(N,N);
            routingMatrixTemp = zeros(N*N,M);
            currentPath = currentDest(k,:);
            currentPath = currentPath(currentPath ~= 0);
            for ell = 2:length(currentPath) % No cycles allowed.
%                 flowsUpdate(currentPath(ell-1),currentPath(ell)) = 1/a;
                routingMatrixTemp((j-1)*N+i,links==currentPath(ell-1)...
                    +(currentPath(ell)-1)*N) = 1/a;
%                 SPlinkNumber(links==currentPath(ell-1)+(currentPath(ell)-1)*N,...
%                     length(currentPath)-1) = SPlinkNumber(links==...
%                     currentPath(ell-1)+(currentPath(ell)-1)*N,...
%                     length(currentPath)-1) + 1;
%                 numberOfSP(links==currentPath(ell-1)+(currentPath(ell)-1)*N,...
%                     cost(k)) = numberOfSP(links==...
%                     currentPath(ell-1)+(currentPath(ell)-1)*N,...
%                     cost(k)) + 1;
%                 numberOfSPW(links==currentPath(ell-1)+(currentPath(ell)-1)*N,...
%                     cost(k)) = numberOfSPW(links==...
%                     currentPath(ell-1)+(currentPath(ell)-1)*N,...
%                     cost(k)) + demand((j-1)*N+i)/a;
            end
            routingMatrix = routingMatrix + routingMatrixTemp;
%             flows = flows + flowsUpdate;
        end
    end
end
% meanNumberOfSP = numberOfSP*[1:N*round(max(max(pathCosts)))]'./sum(numberOfSP,2);
% overallMeanNumberOfSP = sum(numberOfSP,1)*[1:N*round(max(max(pathCosts)))]'/sum(sum(numberOfSP));
% meanNumberOfSPW = numberOfSPW*[1:N*round(max(max(pathCosts)))]'./sum(numberOfSPW,2);
% overallMeanNumberOfSPW = sum(numberOfSPW,1)*[1:N*round(max(max(pathCosts)))]'/sum(sum(numberOfSPW));
% meanSPlinkNumber = SPlinkNumber*[1:N*round(max(max(pathCosts)))]'./sum(SPlinkNumber,2);
% overallMeanSPlinkNumber = sum(SPlinkNumber,1)*[1:N*round(max(max(pathCosts)))]'/sum(sum(SPlinkNumber));

meanNumberOfSP = 1;
overallMeanNumberOfSP = 1;
meanNumberOfSPW = 1;
overallMeanNumberOfSPW = 1;
meanSPlinkNumber = 1;
overallMeanSPlinkNumber = 1;


end