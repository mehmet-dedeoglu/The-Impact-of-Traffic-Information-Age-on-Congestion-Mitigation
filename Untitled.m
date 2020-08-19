
cong1 = zeros(120,3);
tt = 1:120;
reconfig1 = 1;

i = 1;
cong1(:,i) = average_Dep_Net_Congestion1;
update11 = [cong1(1,i);cong1(reconfig1(2:end,1)==1,1)];
ttt1 = [tt(1),tt(reconfig1(2:end,1)==1)];

i = 2;
cong1(:,i) = average_Dep_Net_Congestion1;
update22 = [cong1(1,i);cong1(reconfig1(2:end,1)==1,2)];
ttt2 = [tt(1),tt(reconfig1(2:end,1)==1)];

i = 3;
cong1(:,i) = average_Dep_Net_Congestion1;
update33 = [cong1(1,i);cong1(reconfig1(2:end,1)==1,3)];
ttt3 = [tt(1),tt(reconfig1(2:end,1)==1)];

name = ['cong1 ','update11 ','update22 ','update33 ','ttt1 ','ttt2 ','ttt3'];
save(name,'cong1','update11','update22','update33','ttt1','ttt2','ttt3');

clear all;

figure;plot(tt,cong1(:,1),'b','LineWidth',1.5);hold on;
plot(tt,cong1(:,2),'r','LineWidth',1.5);hold on;
plot(tt,cong1(:,3),'g','LineWidth',1.5);hold on;
plot(ttt1, update11,'b*','MarkerSize',10,'LineWidth',1.5);hold on;
plot(ttt2, update22,'r*','MarkerSize',10,'LineWidth',1.5);hold on;
plot(ttt3, update33,'g*','MarkerSize',10,'LineWidth',1.5);hold on;
xlabel('Iteration');ylabel('Maximum Link Load');grid on;


