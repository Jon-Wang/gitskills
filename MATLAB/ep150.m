%��3.3.2 Kalman�˲��㷨�������
clear;clc;Bushu=50;
Q=0.81; %w(t)�ķ���
R=1;    %V(t)�ķ���
randn('seed',1);w=sqrt(Q)*randn(1,Bushu);
randn('seed',2);v=sqrt(R)*randn(1,Bushu+1);
fai=[0.9 0;-0.6 0.4];gama=[1 2]';H=[1 1];
x(:,1)=[0,0]';y(1)=H*x(:,1)+v(1);
for i=2:Bushu+1
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
%------------------Kalman�˲���----------------%
n=2; %ά��
xjian(:,1)=zeros(2,1);p(:,:)=zeros(2);
for i=1:Bushu
   xxjian(:,i+1)=fai*xjian(:,i);%xxjianΪx��Ԥ��
   e(:,i+1)=y(i+1)-H*xxjian(:,i+1);%��Ϣ
   pp(:,n*(i-1)+1:n*i)=fai*p(:,n*(i-1)+1:n*i)*fai'+gama*Q*gama';
   kf(:,i+1)=pp(:,n*(i-1)+1:n*i)*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'*R);
   xjian(:,i+1)=xxjian(:,i+1)+kf(:,i+1)*e(:,i+1);
   p(:,n*i+1:n*(i+1))=[eye(2)-kf(:,i+1)*H]*pp(:,n*(i-1)+1:n*i);
end
%-------------------��ͼ����--------------------%
t=1:Bushu;
subplot(2,2,1);plot(t,x(1,t),'b',t,xjian(1,t),'r:');
subplot(2,2,2);plot(t,x(2,t),'b',t,xjian(2,t),'r:');


