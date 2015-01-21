%例3.7.3 稳态Kalman滤波算法仿真程序
clc;clear;Bushu=200;
%------------------生成系统模型----------------%
Q=1;R=1;fai=[1.3 1;-0.4 0];gama=[1 -0.5]';H=[1,0];
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
    kf(:,i-1)=temp*H'*inv(H*temp*H'+R);
    kp(:,i-1)=fai*kf(:,i-1);
end
 %-------------------滤波算法--------------------%
 xjian(:,1)=zeros(n,1);pusaif=(eye(n)-kf(:,Bushu)*H)*fai;
 for i=2:Bushu
     xjian(:,i)=pusaif*xjian(:,i-1)+kf(:,Bushu)*y(i);
 end
 %-------------------作图部分--------------------%
 t=1:15;
 subplot(2,2,1);
 plot(t,pp(1,n*(t-1)+1),'r:',t,pp(1,n*(t-1)+2),'k:',t,pp(2,n*(t-1)+2),'b:');
 pp(1,n*(15-1)+1)
 pp(1,n*(15-1)+2)
 pp(2,n*(15-1)+2)

 line([0,15],[1.369,1.369]);line([0,15],[-0.685,-0.685]);
 line([0,15],[0.342,0.342]);%axis([1 15 -1 3]);
 %K图
 subplot(2,2,2);
 plot(t,kf(1,t),'r:',t,kf(2,t),'k:');%axis([1 15 -0.5 1]);
 line([0,15],[0.57805,0.57805]);line([0,15],[-0.28903,-0.28903]);
 %x图
 t=1:200;
 subplot(2,2,3);plot(t,x(1,t),'r:',t,xjian(1,t),'b');axis([0 200 -5 5]);
 subplot(2,2,4);plot(t,x(2,t),'r:',t,xjian(2,t),'b');axis([0 200 -2 2]);
    