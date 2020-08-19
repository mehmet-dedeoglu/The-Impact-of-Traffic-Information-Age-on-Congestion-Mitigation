function reordered=reorder(T,K)

% Inputs
% 
% Variables
% 
% Outputs
% 

j=0;
k=0;
reordered = zeros(2,2*K);
for i=1:2*K
    if min(T(i-j,1),T(K+i-k,1))==T(i-j,1)
        reordered(i,1)=T(i-j,1);
        reordered(i,2)=0;
        k=k+1;
        a=a+[zeros(1,i-1),ones(1,2*K-i+1)];
    else
        reordered(i,1)=T(K+i-k,1);
        reordered(i,2)=1;
        j=j+1;
        a=a-[zeros(1,i-1),ones(1,2*K-i+1)];
    end
end

end