clc;clear;Bushu=4000;
randn('seed',2);e=sqrt(0.25)*randn(1,Bushu);
a1=-1.3;a2=0.4;d1=-0.3;
y(1)=e(1);y(2)=-a1*y(1)+e(2)+d1*e(1);
for i=3:Bushu
        y(i)=-a1*y(i-1)-a2*y(i-2)+e(i)+d1*e(i-1);
end
%赋初值
fai(:,1)=zeros(3,1);fai(:,2)=[-y(1),0,e(1)]';ejian(1)=y(1);
sita(:,1)=zeros(3,1);p(:,1:3)=10^5*eye(3);
for i=2:Bushu
    sita(:,i)=sita(:,i-1)+p(:,3*(i-1)-2:3*(i-1))*fai(:,i)/(1+fai(:,i)'*p(:,3*(i-1)-2:3*(i-1))*fai(:,i))*(y(i)-fai(:,i)'*sita(:,i-1));
    p(:,3*(i)-2:3*(i))=p(:,3*(i-1)-2:3*(i-1))-p(:,3*(i-1)-2:3*(i-1))*fai(:,i)*fai(:,i)'*p(:,3*(i-1)-2:3*(i-1))/[1+fai(:,i)'*p(:,3*(i-1)-2:3*(i-1))*fai(:,i)];
    ejian(i)=y(i)-fai(:,i)'*sita(:,i-1);
    fai(:,i+1)=[-y(i),-y(i-1),ejian(i)]';
end
%噪声方差辨识
taoe(1)=ejian(1)^2;
for i=2:Bushu
    taoe(i)=taoe(i-1)+1/i*(ejian(i)^2-taoe(i-1));
end
figure
t=1:Bushu;%作图部分
subplot(2,2,1);plot(t,sita(1,t),'r');line([0,Bushu],[a1,a1]);
subplot(2,2,2);plot(t,sita(2,t),'r');line([0,Bushu],[a2,a2]);
subplot(2,2,3);plot(t,sita(3,t),'r');line([0,Bushu],[d1,d1]);
subplot(2,2,4);plot(t,taoe(t),'r');line([0,Bushu],[0.25,0.25]);