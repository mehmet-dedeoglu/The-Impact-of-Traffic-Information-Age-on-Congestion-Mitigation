function [count,T_ordered]=Counting_Process(TD,tD,D)

% Inputs
% 
% Variables
% 
% Outputs
% 

K = nnz(TD);
count = zeros(1,2*K);
TD = TD(2:K+1);
TDE = TD+tD;
T = [TD',TDE';zeros(1,K),ones(1,K)];
T_ordered = zeros(2,2*K);

% j=0;
% k=0;
% for i=1:2*K
%     if min(T(1,i-j),T(1,K+i-k))==T(1,i-j)
%         T_ordered(1,i)=T(1,i-j);
%         T_ordered(2,i)=0;
%         k=k+1;
%         count=count+[zeros(1,i-1),ones(1,2*K-i+1)];
%     else
%         T_ordered(1,i)=T(1,K+i-k);
%         T(1,K+i-k) = 2*D;
%         T_ordered(2,i)=1;
%         j=j+1;
%         count=count-[zeros(1,i-1),ones(1,2*K-i+1)];
%     end
% end

% count2 = ones(1,2*K);
% T_ordered2 = zeros(2,2*K);
[T_ordered(1,:),I] = sort(T(1,:),'ascend');
T2 = T(2,:);
T_ordered(2,:) = 2*(1-T2(I))-1;
for i=2:2*K
    count(i) = count(i-1) + T_ordered(2,i);
end


end