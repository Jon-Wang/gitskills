clc;clear;Bushu=2000;lamda=0.99;
randn('seed',4);e=sqrt(0.9^2)*randn(1,Bushu);rand('seed',1)
u=rand(1,Bushu)*2-1;
a(1:Bushu/2)=-0.9;a(Bushu/2+1:Bushu)=0.9;
b(1:Bushu)=3;b(Bushu/2+1:Bushu)=1;
y(1)=e(1);
%y(1)=e(1);y(2)=1.4*y(1)+e(2);
for i=2:Bushu
    y(i)=a(i)*y(i-1)+b(i)*u(i-1)+e(i);
end
%sita(:,Bushu)=zeros(2,1);p(:,1:2)=10^5*eye(2);p(:,3:4)=10^5*eye(2);
fai(:,1)=[0 0]';sita(:,1)=zeros(2,1);
p(:,1:2)=10^5*eye(2);p(:,3:4)=10^5*eye(2);
for i=2:Bushu
    fai(:,i)=[y(i-1),u(i-1)];
    sita(:,i)=sita(:,i-1)+p(:,2*(i-1)-1:2*(i-1))*fai(:,i)/(lamda+fai(:,i)'*p(:,2*(i-1)-1:2*(i-1))...
        *fai(:,i))*(y(i)-fai(:,i)'*sita(:,i-1));
    p(:,2*i-1:2*i)=1/lamda*p(:,2*(i-1)-1:2*(i-1))-p(:,2*(i-1)-1:2*(i-1))*fai(:,i)*fai(:,i)'*p(:,2*(i-1)-1:2*(i-1))/(lamda+fai(:,i)'...
    *p(:,2*(i-1)-1:2*(i-1))*fai(:,i));
end
% for i=1:Bushu
%     emixiu(i)=y(i)-fai(1:2,i)'*sita(1:2,i);
% end
% taoe(1)=emixiu(1)^2;
% for i=2:Bushu
%     taoe(i)=taoe(i-1)+1/i*[emixiu(i)^2-taoe(i-1)];
% end
t=1:Bushu
figure
plot(t,sita(1,t),'r');line([0,Bushu/2],[a(1),a(Bushu/2)]);line([Bushu/2+1,Bushu],[a(Bushu/2+1),a(Bushu)]);
line([Bushu/2,Bushu/2],[a(Bushu/2),a(Bushu/2+1)]);axis([0,Bushu,-2,2]);
figure
plot(t,sita(2,t),'r');line([0,Bushu/2],[b(1),b(Bushu/2)]);line([Bushu/2+1,Bushu],[b(Bushu/2+1),b(Bushu)]);
line([Bushu/2,Bushu/2],[b(Bushu/2),b(Bushu/2+1)]);axis([0,Bushu,0,5]);
