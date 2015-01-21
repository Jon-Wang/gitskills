%习题2.14 有色观测噪声AR(1)模型
clear;clc;Bushu=4000;
randn('seed',1);e=sqrt(1)*randn(1,Bushu);%1,16
randn('seed',16);w=sqrt(0.04)*randn(1,Bushu);
v(1)=w(1); %有色噪声的v(t)生成
for i=2:Bushu
    v(i)=-0.9*v(i-1)+w(i);
end
y(1)=e(1);y(2)=0.8*y(1)+e(2);z(1)=y(1)+v(1);z(2)=y(2)+v(2);
for i=3:Bushu
    y(i)=0.8*y(i-1)-0*y(i-2)+e(i);
    z(i)=y(i)+v(i);
end
Rv0=0.2105263;Rv(1)=-0.9*Rv0;Rv(2)=0.9^2*Rv0;
%Rv0=0.2747252;Rv(1)=0.3*Rv0;Rv(2)=0.3^2*Rv0;
r=[Rv(1),Rv(2)]';
R=[Rv0,Rv(1);Rv(1),Rv0];
Na=2;
fai(:,1)=[0,0]';fai(:,2)=[z(1),0]';
sitab(:,1)=zeros(2,1);sitab(:,2)=zeros(2,1);
p(:,1:2)=10^5*eye(Na);p(:,3:4)=10^5*eye(Na);
sita(:,1)=sitab(:,1);sita(:,2)=sitab(:,2);
for i=Na+1:Bushu
    fai(:,i)=[z(i-1),z(i-2)]';
    sitab(:,i)=sitab(:,i-1)+p(:,Na*(i-2)+1:Na*(i-1))*fai(:,i)/(1+fai(:,i)'*p(:,Na*(i-2)+1:Na*(i-1))*fai(:,i))*(z(i)-fai(:,i)'*sitab(:,i-1));
    p(:,Na*(i-1)+1:Na*i)=p(:,Na*(i-2)+1:Na*(i-1))-p(:,Na*(i-2)+1:Na*(i-1))*fai(:,i)*fai(:,i)'*p(:,Na*(i-2)+1:Na*(i-1))/(1+fai(:,i)'*p(:,Na*(i-2)+1:Na*(i-1))*fai(:,i));
   %--------sita偏差补偿部分-------%
    sita(:,i)=sitab(:,i)-i*p(:,Na*(i-2)+1:Na*(i-1))*r+i*p(:,Na*(i-2)+1:Na*(i-1))*R*sita(:,i-1);
end
%------------估计噪声的方差-----------%
for i=1:Bushu
    emixiu(i)=z(i)-fai(:,i)'*sitab(:,i);
end
taoem(1)=emixiu(1)^2;
for i=2:Bushu
    taoem(i)=taoem(i-1)+1/i*[emixiu(i)^2-taoem(i-1)];
end
%------------taoe偏差补偿部分-----------%
sum=0;
for i=1:Bushu
    for s=1:Na
        for p=s:Na
            temp(s,p)=sita(s,i)*sita(p,i)*R(p-s+1,1);
            sum=sum+temp(s,p);
        end
    end
    taoe(i)=taoem(i)-sum;sum=0;
end
%---------------作图部分----------------%
t=1:Bushu;
figure
subplot(2,2,1);plot(t,sitab(1,t),'r');line([0,Bushu],[0.8,0.8]);axis([0 Bushu -2 2]);
subplot(2,2,2);plot(t,sitab(2,t),'r');line([0,Bushu],[0,0]);axis([0 Bushu -1 1]);
subplot(2,2,3);plot(t,taoem,'r');line([0,Bushu],[1,1]);axis([0 Bushu -2 2]);
figure
subplot(2,2,1);plot(t,sita(1,t),'r');line([0,Bushu],[0.8,0.8]);axis([0 Bushu -2 2]);
subplot(2,2,2);plot(t,sita(2,t),'r');line([0,Bushu],[0,0]);axis([0 Bushu -1 1]);
subplot(2,2,3);plot(t,taoe,'r');line([0,Bushu],[1,1]);axis([0 Bushu -2 2]);



