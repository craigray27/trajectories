c2=ForwardEulerSolver;
c3=SymplecticEulerSolver;
c4=StoermerVerletSolver;
q=zeros(10001,1);
z=zeros(10001,1);
c=zeros(10001,1);
for i=1:10001
    q(i)=(sqrt(c2(i,1)^2+c2(i,2)^2)-3.844e8)^2;
    z(i)=(sqrt(c4(i,1)^2+c4(i,2)^2)-3.844e8)^2;
    c(i)=(sqrt(c3(i,1)^2+c3(i,2)^2)-3.844e8)^2;
end
e1=sum(q)/10001;
e2=sum(z)/10001;
e3=sum(c)/10001;