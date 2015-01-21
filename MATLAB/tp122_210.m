clear all;clc;Bushu=3000;
% * * * * * * 生成模型 * * * * * * %
a1=1.2;a2=-0.35;Qe=1;Qv=0.16;
randn('seed',1);e=sqrt(Qe)*randn(1,Bushu+10);
randn('seed',9);v=sqrt(Qv)*randn(1,Bushu+10);
y(1)=e(1);y(2)=a1*y(1)+e(2);
y(Bushu+10)=0; %扩维，节约运算时间
z(1)=y(1)+v(1);z(2)=y(2)+v(2);z(Bushu+10)=0;
for i=3:Bushu+10
    y(i)=a1*y(i-1)+a2*y(i-2)+e(i);
    z(i)=y(i)+v(i);
end
% * * * * * * 系统参数辨识 * * * * * * %
Na=2; %因为待估参数有两个
fai(:,Bushu)=[0,0]';faijian(:,Bushu)=[0,0]';sita(:,Bushu)=zeros(2,1);
p(:,1:2)=10^5*eye(Na);p(:,3:4)=10^5*eye(Na);
p(:,Na*(Bushu-1)+1:Na*Bushu)=10^5*eye(Na);
for i=Na+1:Bushu
    fai(:,i)=[z(i-1),z(i-2)]';
    faijian(:,i)=fai(:,i-2);
sita(:,i)=sita(:,i-1)+p(:,Na*(i-2)+1:Na*(i-1))*faijian(:,i)/[1+fai(:,i)'*p(:,Na*(i-2)+1:Na*(i	-1))*faijian(:,i)]*(z(i)-fai(:,i)'*sita(:,i-1));
p(:,Na*(i-1)+1:Na*i)=p(:,Na*(i-2)+1:Na*(i-1))-p(:,Na*(i-2)+1:Na*(i-1))*faijian(:,i)/(1+	fai(:,i)'*p(:,Na*(i-2)+1:Na*(i-1))*faijian(:,i))*...
        fai(:,i)'*p(:,Na*(i-2)+1:Na*(i-1));
end
% * * * * * * 噪声方差的估计 * * * * * * %
zyiheng(1)=z(1);zyiheng(2)=z(2)-sita(1,2)*z(1);zyiheng(Bushu)=0;
for i=3:Bushu
    zyiheng(i)=z(i)-sita(1,i)*z(i-1)-sita(2,i)*z(i-2);
end
R0=0;R(2,1)=R0+1/1*[0-R0];R(2,2)=R(2,1)+1/2*[0-R(2,1)];
taov(1)=-R(2,1)/sita(2,1);tao(2)=-R(2,2)/sita(2,2);
for t=3:Bushu
    R(2,t)=R(2,t-1)+1/t*[zyiheng(t-2)*zyiheng(t)-R(2,t-1)];
    taov(t)=-R(2,t)/sita(2,t);
end
taoz0(Bushu)=0;
for i=2:Bushu
    taoz0(i)=taoz0(i-1)+1/i*[zyiheng(i)^2-taoz0(i-1)];
end
for i=1:Bushu
    aj(i)=1+sita(1,i)^2+sita(2,i)^2;
end
for i=1:Bushu
    taoem(i)=taoz0(i)-taov(1,t)*aj(i);
end
% * * * * * * 作图部分 * * * * * * %
t=1:Bushu;
subplot(2,2,1);plot(sita(1,:),'r');line([0,Bushu],[a1,a1]);axis([0 Bushu -1 4]);xlabel('t/step');
subplot(2,2,2);plot(sita(2,:),'r');line([0,Bushu],[a2,a2]);axis([0 Bushu -2 2]);
subplot(2,2,3);plot(taov,'r');line([0,Bushu],[Qv,Qv]);axis([0 Bushu -2 6]);
subplot(2,2,4);plot(taoem,'r');line([0,Bushu],[Qe,Qe]);axis([0 Bushu -5 20]);
sita(1,Bushu)
sita(2,Bushu)
taov(Bushu)
taoem(Bushu)

