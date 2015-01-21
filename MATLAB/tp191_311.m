%习题3.11 稳态Kalman滤波算法仿真程序
%path(path,'D:\shortcut\MATLAB')
clear all;clc;Bushu=300;
fai=[1 1;0 0];gama=[1 1]';H=[1,0];R=0.04;
n=2;%维数

Taog=9;Lmd=0.1;Q=Lmd*Taog;
randn('seed',5);g=sqrt(Taog)*randn(1,Bushu+10);
rand('state',1);k=rand(1,Bushu+10);
for i=1:Bushu+10;
    if k(i)<=Lmd;b(i)=1;else;b(i)=0;end;
    w(i)=b(i)*g(i);
end

randn('seed',3);v=sqrt(R)*randn(1,Bushu+10);
x(:,1)=zeros(2,1);y(1,1)=H*x(:,1)+v(1);
for i=2:Bushu+10;
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
%----------------------------------%
[xigema,kf,kp,pusaif,pusaip,e,xjian,xxjian]=TS_WTKalman(fai,gama,H,Q,R,y);
N=3;
Qe=H*xigema*H'+R;
 for i=1:N
     M(i,1)=Q*gama'*pusaip'^(i-1)*H'*inv(Qe);
 end
  for i=1:Bushu
     wjian(1,i)=M(1,1)*e(i+1);
     wjian(2,i)=wjian(1,i)+M(2,1)*e(i+2);
     wjian(3,i)=wjian(2,i)+M(3,1)*e(i+3);
 end
 %-------------------作图部分--------------------%
 figure
 t=1:Bushu;
 subplot(1,2,1);plot(t,wjian(1,t),'k.');%axis([0 300 -2.5 2.5])
 for t=1:Bushu;hh=line([t,t],[0,w(t)]);set(hh,'color','k');end
 t=1:Bushu;
 subplot(1,2,2);plot(t,wjian(2,t),'k.');%axis([0 300 -2.5 2.5])
 for t=1:Bushu;hh=line([t,t],[0,w(t)]);set(hh,'color','k');end
 figure
 t=1:Bushu;
 subplot(1,1,1);plot(t,wjian(3,t),'k.');%axis([0 300 -2.5 2.5])
 for t=1:Bushu;hh=line([t,t],[0,w(t)]);set(hh,'color','k');end
