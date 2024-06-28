
plot3(bodies(:,1),bodies(:,2),bodies(:,3))
hold on
plot3(bodies(:,4),bodies(:,5),bodies(:,6))
hold on
plot3(bodies(:,7),bodies(:,8),bodies(:,9))
hold on
plot3(bodies(:,10),bodies(:,11),bodies(:,12))
hold on
plot3(bodies(:,13),bodies(:,14),bodies(:,15))
xlabel('x')
ylabel('y')
zlabel('z')
legend('Sun','Mercury','Venus','Earth','Mars')
title('trajectories of five bodies')
c5=zeros(20,15);
k=1;
for i=1:500:10000
    c5(k,:)=bodies(i,:);
    k=k+1;
end
c5=c5(:,1:9);
