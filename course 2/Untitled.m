Eh=StoermerVerletSolverEh;
d=Eh(:,2:3);
x=d(:,1);
y=log(d(:,2));
hold on;
z=x+0.1187;
plot(d(:,1),z,'r');
E=zeros(30,3);
k=1;
for i=1:33:990
E(k,:)=Eh(i,:);
k=k+1;
end

