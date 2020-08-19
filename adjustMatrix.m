function [P,pindnew] = adjustMatrix(P,leng,i)

pindnew = zeros(size(i));
for j = 1:length(i)
    pindnew(j) = i(j)+(leng-1)*(j-1);
    Ptemp = P(pindnew(j),:).*ones(leng,1);
    P = [P(1:pindnew(j)-1,:);Ptemp;P(pindnew(j)+1:end,:)];
end

end