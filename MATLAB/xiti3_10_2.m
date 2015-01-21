clc;clear;Bushu=200;
% * * * * * * 生成系统模型 * * * * * * %
Q=1;R=1;fai=0.5;gama=1;H=1;
randn('seed',4);w=sqrt(Q)*randn(1,Bushu+10);
randn('seed',1);v=sqrt(R)*randn(1,Bushu+10);
x(:,1)=0;
y(1)=H*x(:,1)+v(1);
for i=2:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
% * * * * * * 递推Riccati方程 * * * * * * %
n=1; %维数
pp(:,:)=10^5*eye(n);
for i=2:Bushu+1;
    temp=pp(:,n*(i-2)+1:n*(i-1));
    kf(:,i-1)=temp*H'*inv(H*temp*H'+R);
    pp(:,n*(i-1)+1:n*i)=fai*[temp-temp*H'*inv(H*temp*H'+R)*H*temp]*fai'+gama*Q*gama';
    kp(:,i-1)=fai*kf(:,i-1);
end
 % * * * * * * 稳态滤波和预报 * * * * * * %
 xjian(:,1)=zeros(n,1);xxjian(:,1)=zeros(n,1);e(1)=y(1)-H*xxjian(:,1);
 pusaif=(eye(n)-kf(:,Bushu)*H)*fai;
 pusaip=fai*(eye(n)-kf(:,Bushu)*H);
 Qe=H*temp*H'+R;
 for i=2:Bushu+10
     xjian(:,i)=pusaif*xjian(:,i-1)+kf(:,Bushu)*y(i);
     e(i-1)=y(i-1)-H*xxjian(:,i-1);
     xxjian(:,i)=fai*xxjian(:,i-1)+kp(:,Bushu)*e(i-1);
     xxnjian(:,i)=fai*xxjian(:,i); %2步预报
 end
 % * * * * * * 输入白噪声估值器 * * * * * * %
for i=1:Bushu+5
    M(1,i)=Q*gama'*H'*inv(Qe);
    M(2,i)=Q*gama'*pusaip'*H'*inv(Qe);
end
for i=1:Bushu
    wjian(1,i)=M(1,i+1)*e(i+1);             %一步平滑
    wjian(2,i)=wjian(1,i)+M(2,i+2)*e(i+2);  %二步平滑
end
% * * * * * * 观测白噪声估值器 * * * * * * %
for i=1:Bushu+5
    Mv(1,i)=R*inv(Qe);
    Mv(2,i)=-R*kf(:,Bushu)'*fai'*H'*inv(Qe);
    Mv(3,i)=-R*kf(:,Bushu)'*fai'*pusaip'*H'*inv(Qe);
end
for i=1:Bushu
    vjian(1,i)=Mv(1,i)*e(i);                 %一步平滑
    vjian(2,i)=vjian(1,i)+Mv(2,i+1)*e(i+1);  %二步平滑
    vjian(3,i)=vjian(2,i)+Mv(3,i+2)*e(i+2);  %三步平滑
end
 figure
 t=1:Bushu;
 subplot(2,2,1);plot(t,x(1,t),t,xjian(1,t),'r:');  %稳态Kalman滤波器
 subplot(2,2,2);plot(t,x(1,t),t,xxnjian(1,t),'r:');%稳态Kalman**2步**预报器
 subplot(2,2,3);plot(t,w(1,t),t,wjian(1,t),'r:');  %输入白噪声1步预报
 subplot(2,2,4);plot(t,w(1,t),t,wjian(2,t),'r:');  %输入白噪声2步预报
 figure
 t=1:Bushu;
 subplot(1,2,1);plot(t,v(1,t),t,vjian(1,t),'r:');  %观测白噪声0步预报
 subplot(1,2,2);plot(t,v(1,t),t,vjian(2,t),'r:');  %观测白噪声1步预报
  figure
 subplot(1,1,1);plot(t,v(1,t),t,vjian(3,t),'r:');  %观测白噪声2步预报
 