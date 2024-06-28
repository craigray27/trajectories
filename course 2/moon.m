G=6.674e-11;
masse=5.972e24;
massm=7.342e22;
r=3.844e8;
T=2*pi*sqrt(r^3/(G*masse));
v0=sqrt(G*masse/r);
ve=[0,0,0];
vm=[0,v0,0];
h=T/10000;
xe=[0,0,0];
xm=[r,0,0];
xearth=zeros(10000,3);
xmoon=zeros(10000,3);
d=zeros(10000,1);
fee=zeros(10000,3);
fmm=zeros(10000,3);
for i=1:10000
xearth(i,:)=xe;
xmoon(i,:)=xm;
fe=(xm-xe)*(G*massm/(norm(xm-xe))^3);
fee(i,:)=fe;
newxe=xe+h*ve+fe*h^2/2;
newfe=(xm-newxe)*(G*massm/(norm(xm-newxe))^3);
newve=ve+(fe+newfe)*h/2;
ve=newve;
xe=newxe;

fm=(xe-xm)*(G*masse/(norm(xm-xe))^3);
fmm(i,:)=fm;
newxm=xm+h*vm+fm*h^2/2;
newfm=(xe-newxm)*(G*masse/(norm(xm-newxe))^3);
newvm=vm+(fm+newfm)*h/2;
vm=newvm;
xm=newxm;

d(i)=norm(xe-xm);

end





