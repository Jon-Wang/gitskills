%习题3.9 稳态Kalman滤波算法仿真程序
clc;clear;Bushu=200;
%------------------生成系统模型----------------%
Q=1;R=1;fai=[1 1;0 0];gama=[1 1]';H=[1,0];
randn('seed',1);w=sqrt(Q)*randn(2,Bushu+10);
randn('seed',2);v=sqrt(R)*randn(1,Bushu+10);
x(:,1)=[0 0]';
y(1)=H*x(:,1)+v(1);
for i=2:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
%------------------递推Riccati方程----------------%
n=2;%维数
pp(:,:)=eye(n);
for i=2:Bushu+1;
    temp=pp(:,n*(i-2)+1:n*(i-1));
    pp(:,n*(i-1)+1:n*i)=fai*[temp-temp*H'*inv(H*temp*H'+R)*H*temp]*fai'+gama*Q*gama';
    kf(:,i-1:i-1)=temp*H'*inv(H*temp*H'+R);
    kp(:,i-1:i-1)=fai*kf(:,i-1:i-1);
end
 %-------------------滤波算法--------------------%
 xjian(:,1)=zeros(n,1);xxjian(:,1)=zeros(n,1);
 pusaif=(eye(n)-kf(:,Bushu:Bushu)*H)*fai;
 pusaip=fai*(eye(n)-kf(:,Bushu:Bushu)*H);
 for i=2:Bushu+1
     %xxjian(:,i+1)=fai*xjian(:,i);%xxjian为x的预报
     xjian(:,i)=pusaif*xjian(:,i-1)+kf(:,Bushu:Bushu)*y(i);
     xxjian(:,i)=pusaip*xxjian(:,i-1)+kp(:,Bushu:Bushu)*y(i);
 end
 %-------------------作图部分--------------------%
 figure
 t=1:200;
 subplot(2,2,1);plot(t,x(1,t),'r:',t,xjian(1,t),'b');%axis([0 200 -5 5]);
 subplot(2,2,2);plot(t,x(2,t),'r:',t,xjian(2,t),'b');%axis([0 200 -2 2]);
 subplot(2,2,3);plot(t,x(1,t),'r:',t,xxjian(1,t),'b');%axis([0 200 -5 5]);
 subplot(2,2,4);plot(t,x(2,t),'r:',t,xxjian(2,t),'b');%axis([0 200 -2 2]);

