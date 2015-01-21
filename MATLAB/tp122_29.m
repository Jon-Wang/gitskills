clc;clear;Bushu=3000;
Q=0.36;a1=-1.5;a2=0.54;d1=-1.2;d2=0.35;
randn('seed',1);e=sqrt(Q)*randn(1,Bushu);
%生成模型
y(1)=e(1);y(2)=-a1*y(1)+e(2)+d1*e(1);y(3)=-a1*y(2)-a2*y(1)+e(3)+d1*e(2)+d2*e(1);
for i=4:Bushu
        y(i)=-a1*y(i-1)-a2*y(i-2)+e(i)+d1*e(i-1)+d2*e(i-2);
end
%拟合高阶AR(15)
n=15;
fai(:,1)=zeros(1,n)';
for i=2:n
    for j=1:i-1
        fai(j,i)=-y(i-j);
    end
end
sita(:,1)=zeros(1,n)';p(:,1:n)=10^5*eye(n);
for i=2:Bushu
    if i>15
        for j=1:n
            fai(j,i)=-y(i-j);
        end
    end        
    sita(:,i)=sita(:,i-1)+p(:,n*(i-2)+1:n*(i-1))*fai(:,i)/[1+fai(:,i)'*p(:,n*(i-2)+1:n*(i-1))*fai(:,i)]*[y(i)-fai(:,i)'*sita(:,i-1)];
    p(:,n*(i-1)+1:n*i)=p(:,n*(i-2)+1:n*(i-1))-p(:,n*(i-2)+1:n*(i-1))*fai(:,i)*fai(:,i)'*p(:,n*(i-2)+1:n*(i-1))/[1+fai(:,i)'*p(:,n*(i-2)+1:n*(i-1))*fai(:,i)];
end
%下面用LS算法
for i=1:Bushu
    B(:,1)=[sita(:,i);0];
    B(:,2)=[1;sita(:,i)];
    b=[-sita(2:15,i);0;0];
    djian(:,i)=inv(B'*B)*B'*b;
    dheng=[djian(1,i)]';
    D=[1];
    bheng=[sita(1,i)]';
    ajian(:,i)=dheng+D*bheng;
end
%下面估计噪声的方差
for i=1:Bushu
    emixiu(i)=y(i)-fai(:,i)'*sita(:,i);
end
taoe(1)=emixiu(1)^2;
for i=2:Bushu
    taoe(i)=taoe(i-1)+1/i*[emixiu(i)^2-taoe(i-1)];
end
%作图部分
t=1:Bushu;
subplot(2,2,1);plot(t,ajian(1,t),'r');line([0,Bushu],[a1,a1]);
subplot(2,2,2);plot(t,ajian(2,t),'r');line([0,Bushu],[a2,a2]);
subplot(2,2,3);plot(t,djian(1,t),'r');line([0,Bushu],[d1,d1]);
subplot(2,2,4);plot(t,djian(2,t),'r');line([0,Bushu],[d2,d2]);
figure;
t=1:Bushu;plot(t,taoe(t),'r');line([0,Bushu],[Q,Q]);


