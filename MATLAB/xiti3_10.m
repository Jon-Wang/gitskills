clc;clear;Bushu=200;
%**************生成系统模型***************%
Q=1;R=1;fai=0.5;gama=1;H=1;
randn('seed',4);w=sqrt(Q)*randn(1,Bushu+10);
randn('seed',1);v=sqrt(R)*randn(1,Bushu+10);
x(:,1)=0;
y(1)=H*x(:,1)+v(1);
for i=2:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
%**************递推Riccati方程***************%
n=1;%维数
pp(:,:)=10^5*eye(n);
for i=2:Bushu+1;
    temp=pp(:,n*(i-2)+1:n*(i-1));
    kf(:,i-1)=temp*H'*inv(H*temp*H'+R);
    pp(:,n*(i-1)+1:n*i)=fai*[temp-temp*H'*inv(H*temp*H'+R)*H*temp]*fai'+gama*Q*gama';
    kp(:,i-1)=fai*kf(:,i-1);
end
 %***********滤波和预报**************%
 xjian(:,1)=zeros(n,1);xxjian(:,1)=zeros(n,1);e(1)=y(1)-H*xxjian(:,1);
 pusaif=(eye(n)-kf(:,Bushu)*H)*fai;pusaip=fai*(eye(n)-kf(:,Bushu)*H);
 for i=2:Bushu
     xjian(:,i)=pusaif*xjian(:,i-1)+kf(:,Bushu)*y(i);
     %xxjian(:,i)=pusaip*xxjian(:,i-1)+kp(:,Bushu)*y(i);
     e(i-1)=y(i-1)-H*xxjian(:,i-1);
     xxjian(:,i)=fai*xxjian(:,i-1)+kp(:,Bushu)*e(i-1);
     xxnjian(:,i)=fai*xxjian(:,i);
 end
 figure
 t=1:Bushu;
 subplot(2,2,1);plot(t,x(1,t),t,xjian(1,t),'r:');
 subplot(2,2,2);plot(t,x(1,t),t,xxnjian(1,t),'r:');
%  subplot(2,2,3);plot(t,x(1,t),t,xxjian(1,t),'r:');
%  subplot(2,2,4);plot(t,x(2,t),t,xxjian(2,t),'r:');