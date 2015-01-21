%习题3.7 稳态Kalman滤波算法仿真程序
clc;clear;Bushu=100;
%------------------生成系统模型----------------%
Q=0.81;R=1;fai=[0.9 0;-0.5 0.2];gama=[1 2]';H=[1,1];
randn('seed',1);w=sqrt(Q)*randn(2,Bushu+10);
randn('seed',1);v=sqrt(R)*randn(1,Bushu+10);
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
    kf(:,i-1)=temp*H'*inv(H*temp*H'+R);
    kp(:,i-1)=fai*kf(:,i-1);
end
 %-------------------滤波算法--------------------%
 xjian(:,1)=zeros(n,1);pusaif=(eye(n)-kf(:,Bushu)*H)*fai;
 for i=2:Bushu
     xjian(:,i)=pusaif*xjian(:,i-1)+kf(:,Bushu)*y(i);
 end
 %-------------------作图部分--------------------%
 t=1:100;
 subplot(2,2,1);
 plot(t,pp(1,n*(t-1)+1),'r:',t,pp(1,n*(t-1)+2),'k:',t,pp(2,n*(t-1)+2),'b:');
 pp(1,n*(100-1)+1)
 pp(1,n*(100-1)+2)
 pp(2,n*(100-1)+2)
 line([0,30],[1.0046,1.0046]);line([0,30],[1.5230,1.5230]);
 line([0,30],[3.3083,3.3083]);axis([1 30 0 4]);
 %K图
 subplot(2,2,2);
 plot(t,kf(1,t),'r:',t,kf(2,t),'k:');axis([1 30 0 0.8]);
 line([0,30],[0.3024,0.3024]);line([0,30],[0.5780,0.5780]);
 %x图
 t=1:Bushu;
 subplot(2,2,3);plot(t,x(1,t),'r:',t,xjian(1,t),'b');%axis([0 Bushu -5 5]);
 subplot(2,2,4);plot(t,x(2,t),'r:',t,xjian(2,t),'b');%axis([0 Bushu -2 2]);
    