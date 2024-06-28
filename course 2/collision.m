G=6.674e-11;
m=5.972e24;
xp=[0,0,0];
r=3.844e8;
r2=6.378e6+1.737e6;
T=2*pi*sqrt(r^3/(G*m));
v0=sqrt(G*m/r);
v=[0,0,0];
h=T/10000;
x=[r,0,0];
t=0;
while(norm(xp-x)>r2)
f=(xp-x)*(G*m/(norm(xp-x))^3);
x=x+v*h;
v=v+f*h;
t=t+h;
end

tf=t;


t=t-h;
v=v-f*h;
x=x-v*h;
h=0.1;
while(norm(xp-x)>r2)
f=(xp-x)*(G*m/(norm(xp-x))^3);
x=x+v*h;
v=v+f*h;
t=t+h;
end
%������16a:0.65̫��ֱ����0.69̫��������4450�档
%������16b:0.23��̫��ֱ����0.203��̫��������3218�档
%������16c:0.75��ľ��ֱ����0.333��ľ���������볤��0.7au���������228�졣