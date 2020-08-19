function [paths,pathCost] = dijkstra(graph_index,weight_index,source,destination,N)

s = source;
d = destination;
wi = weight_index;
gi = graph_index;
[m,~] = size(gi);
m = 10*m;

visited = zeros(1,m);
P = zeros(1,m);
P(:,1) = s;
J = Inf(m);
J(s,:) = zeros(1,m);
k = s;
ell = 2;
while ell ~= 0
    P(:,ell) = P(:,ell-1);
    J(:,ell) = J(:,ell-1);
    pind = find(P(:,ell) == k);
    i = 1;
    leng = length(find(gi(k,:)<N+1));
    [P,pindnew] = adjustMatrix(P,leng,pind);
    for j = 1:leng
        if visited(gi(k,i)) == 0
            if wi(k,i)+J(k,ell-1) < J(gi(k,i),ell)
                deleteInd = find(P(:,ell) == gi(k,i));
                P(deleteInd,:) = [];
                pindnewt = pindnew;
                for cc = 1:length(pindnew)
                    for dd = 1:length(deleteInd)
                        if pindnew(cc) > deleteInd(dd)
                            pindnewt(cc) = pindnewt(cc) - 1;
                        end
                    end
                end
                pindnew = pindnewt;
%                 prevDelLeng = ((deleteInd < (pindnew)));%!!!!!!!
%                 pindnew = pindnew-prevDelLeng;
%                 prevDelLeng = length(deleteInd(deleteInd < min(pindnew)));
%                 pindnew = pindnew - prevDelLeng;
                J(gi(k,i),ell) = wi(k,i)+J(k,ell-1);
                P(pindnew,ell) = gi(k,i);
                pindnew = pindnew + 1;
            elseif wi(k,i)+J(k,ell-1) == J(gi(k,i),ell)
                J(gi(k,i),ell) = wi(k,i)+J(k,ell-1);
                P(pindnew,ell) = gi(k,i);
                pindnew = pindnew + 1;
            else
                P(pindnew,:) = [];
                pindnew = pindnew-((1:length(pindnew))-1)';
            end
        else
            P(pindnew,:) = [];
            pindnew = pindnew-((1:length(pindnew))-1)';
        end
        i = i+1;
    end
    visited(k) = 1;
    next = J(:,ell).*(ones(1,m)-visited)';
    next(next == 0) = inf;
    [~,k]= min(next);
    ell = ell + 1;
    if k == d
        pathInd = find(P(:,ell-1) == d);
        ell = 0;
        pathsTemp = P(pathInd,:);
        paths = zeros(1,N);
        [aa,bb] = size(pathsTemp);
        for i = 1:aa
            path = zeros(1,bb);
            path(1) = pathsTemp(i,1);
            k = 2;
            for j = 2:bb
                if pathsTemp(i,j) ~= pathsTemp(i,j-1)
                    path(k) = pathsTemp(i,j);
                    k = k+1;
                end
            end
            inv = path(path ~= 0);
            asd = inv(end:-1:1);
            paths(i,:) = [asd,path(length(asd)+1:N)];
        end
    end
end
pathCost = next(d);
end