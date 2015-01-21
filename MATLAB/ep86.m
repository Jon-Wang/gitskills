clc;clear;Bushu=1000;
randn('seed',2);e=sqrt(1)*randn(1,Bushu);
randn('seed',6);u=sqrt(1)*randn(1,Bushu);
%生成模型
a1=-0.9;a2=0.2;b2=0.5;d1=-0.4;
y(1)=e(1);y(2)=-a1*y(1)+e(2)+d1*e(1);
for i=3:Bushu
        y(i)=-a1*y(i-1)-a2*y(i-2)+b2*u(i-2)+e(i)+d1*e(i-1);
end
%RLS算法，其目的是得到对噪声的估计
n=15;%拟合高阶AR(15)
fai(:,1)=zeros(1,2*n)';
for i=2:n
    for j=1:i-1
        fai(j,i)=-y(i-j);
        fai(j+n,i)=u(i-j);
    end
end
sita2(:,1)=zeros(1,2*n)';p(:,1:2*n)=10^5*eye(2*n);
for i=2:Bushu
    if i>15
        for j=1:n
            fai(j,i)=-y(i-j);
            fai(j+n,i)=u(i-j);
        end
    end        
    sita2(:,i)=sita2(:,i-1)+p(:,2*n*(i-2)+1:2*n*(i-1))*fai(:,i)/[1+fai(:,i)'*p(:,2*n*(i-2)+1:2*n*(i-1))*fai(:,i)]*[y(i)-fai(:,i)'*sita2(:,i-1)];
    p(:,2*n*(i-1)+1:2*n*i)=p(:,2*n*(i-2)+1:2*n*(i-1))-p(:,2*n*(i-2)+1:2*n*(i-1))*fai(:,i)*fai(:,i)'*p(:,2*n*(i-2)+1:2*n*(i-1))/[1+fai(:,i)'*p(:,2*n*(i-2)+1:2*n*(i-1))*fai(:,i)];
    ejian(i)=y(i)-fai(:,i)'*sita2(:,i);
end
%下面估计噪声的方差
taoe(1)=ejian(1)^2;
for i=2:Bushu
    taoe(i)=taoe(i-1)+1/i*(ejian(i)^2-taoe(i-1));
end
%第一阶拟合结束
%RELS算法
fai2(1:5,1)=zeros(1,5)';fai2(1:5,2)=[y(1),0,u(1),0,ejian(1)]';
sita(1:5,1)=zeros(1,5)';sita(:,2)=zeros(1,5)';
p2(1:5,1:5)=10^5*eye(5);p2(:,6:10)=10^5*eye(5);
for i=3:Bushu
    fai2(:,i)=[-y(i-1),-y(i-2),u(i-1),u(i-2),ejian(i-1)]';
    sita(:,i)=sita(:,i-1)+p2(:,5*(i-2)+1:5*(i-1))*fai2(:,i)*[y(i)-fai2(:,i)'*sita(:,i-1)]/[1+fai2(:,i)'*p2(:,5*(i-2)+1:5*(i-1))*fai2(:,i)];
    p2(:,5*(i-1)+1:5*i)=p2(:,5*(i-2)+1:5*(i-1))-p2(:,5*(i-2)+1:5*(i-1))*fai2(:,i)*fai2(:,i)'*p2(:,5*(i-2)+1:5*(i-1))/[1+fai2(:,i)'*p2(:,5*(i-2)+1:5*(i-1))*fai2(:,i)];
end
%作图部分
t=1:Bushu;
subplot(2,2,1);plot(t,sita(1,t),'r');line([0,Bushu],[a1,a1]);axis([0,Bushu,-1.5,0]);
subplot(2,2,2);plot(t,sita(2,t),'r');line([0,Bushu],[a2,a2]);axis([0,Bushu,-1,1]);
subplot(2,2,3);plot(t,sita(4,t),'r');line([0,Bushu],[b2,b2]);axis([0,Bushu,-1,1]);
subplot(2,2,4);plot(t,sita(5,t),'r');line([0,Bushu],[d1,d1]);axis([0,Bushu,-1,1]);
figure
subplot(2,2,1);plot(t,taoe(t),'r');line([0,Bushu],[1,1]);axis([0,Bushu,0,2]);
