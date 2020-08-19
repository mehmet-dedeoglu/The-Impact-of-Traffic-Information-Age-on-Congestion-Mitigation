function [TD,TO] = Constant_Traffic(N,iN,dT)

% Generate a constant traffic migrating nodes at certain time instants
time_schedule = [0,round(iN/4),2*round(iN/4),3*round(iN/4),iN];
% traffic_amount1 = 5*[1,0,1,0];
% traffic_amount2 = 5*[0,1,0,1];
traffic_amount1 = 5*[1,0,1,0];
traffic_amount2 = 5*[0,1,0,1];


TO = zeros(N,N,iN+2);
TD = zeros(N,N,iN+2);
for k = 1:length(time_schedule)-1
    TD(1:floor(N/2),:,time_schedule(k)+1:time_schedule(k+1)+1) = traffic_amount1(k);
    TD(floor(N/2)+1:end,:,time_schedule(k)+1:time_schedule(k+1)+1) = traffic_amount2(k);
    TO(:,:,time_schedule(k)+1:time_schedule(k+1)+1) = dT*time_schedule(k);
end

end