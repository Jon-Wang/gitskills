clc;clear;Bushu=2000;
randn('seed',1);e=sqrt(0.25)*randn(1,Bushu);
%生成模型
y(1)=e(1);y(2)=e(2)-1.2*e(1);y(3)=e(3)-1.2*e(2)+0.47*e(1);
for i=4:Bushu
        y(i)=e(i)-1.2*e(i-1)+0.47*e(i-2)-0.06*e(i-3);
end
n=15;%拟合高阶AR(15)
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
    for j=1:3
        B(j,j)=1,s=0;
        for k=j+1:15+j
            s=s+1;
            B(k,j)=sita(s,i);
        end
    end
    for k=1:15+3
        if k<=15
            b(k)=-sita(k,i);
        else
            b(k)=0;
        end
    end
    djian(:,i)=inv(B'*B)*B'*b';
end
%下面估计噪声方差
for i=1:Bushu
    emixiu(i)=y(i)-fai(:,i)'*sita(:,i);
end
taoe(1)=emixiu(1)^2;
for i=2:Bushu
    taoe(i)=taoe(i-1)+1/i*[emixiu(i)^2-taoe(i-1)];
end
%作图部分
t=1:Bushu;
subplot(2,2,1);plot(t,djian(1,t),'r');line([0,Bushu],[-1.2,-1.2]);
subplot(2,2,2);plot(t,djian(2,t),'r');line([0,Bushu],[0.47,0.47]);
subplot(2,2,3);plot(t,djian(3,t),'r');line([0,Bushu],[-0.06,-0.06]);
subplot(2,2,4);plot(t,taoe(t),'r');line([0,Bushu],[0.25,0.25]);axis([0,Bushu,0,0.4]);
