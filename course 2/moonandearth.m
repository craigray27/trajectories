G=6.674e-11;
mass=1.9891e30;
r=5.791e10;
T=2*pi*sqrt(r^3/(G*mass));
T2=T/(3600*24);
c4=earthandmoon(:,1:6);
plot3(c4(:,1),c4(:,2),c4(:,3))
hold on
plot3(c4(:,4),c4(:,5),c4(:,6))
xlabel('x')
ylabel('y')
zlabel('z')
legend('earth','moon')
title('trajectories of moon and earth')
c5=zeros(20,6);
k=1;
for i=1:500:10001
    c5(k,:)=c4(i,:);
    k=k+1;
end