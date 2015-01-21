%ϰ��3.3 Kalman�˲��㷨�������
clear;clc;Bushu=300;
Q=0.64; %w(t)�ķ���
R=0.36;    %V(t)�ķ���
randn('seed',1);w=sqrt(Q)*randn(1,Bushu+2);
randn('seed',2);v=sqrt(R)*randn(1,Bushu+2);
fai=[0.6 -0.4;0.2 1.1];gama=[1 2]';H=[1 1];
x(:,1)=[0,0]';y(1)=H*x(:,1)+v(1);
for i=2:Bushu+2
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
%------------------Kalman�˲���----------------%
n=2; %ά��
xjian(:,1)=zeros(2,1);p(:,:)=zeros(2);
for i=1:Bushu+1
   xxjian(:,i+1)=fai*xjian(:,i);%xxjianΪx��Ԥ��
   xxxjian(:,i+1)=fai^2*xjian(:,i);%����Ԥ��
   e(:,i+1)=y(i+1)-H*xxjian(:,i+1);%��Ϣ
   pp(:,n*(i-1)+1:n*i)=fai*p(:,n*(i-1)+1:n*i)*fai'+gama*Q*gama';
   kf(:,i+1)=pp(:,n*(i-1)+1:n*i)*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'+R);
   xjian(:,i+1)=xxjian(:,i+1)+kf(:,i+1)*e(:,i+1);
   p(:,n*i+1:n*(i+1))=[eye(2)-kf(:,i+1)*H]*pp(:,n*(i-1)+1:n*i);
   %mxjian(:,i+1)=xjian(:,i)+pp(:,n*(i-2)+1:n*(i-1))*[eye(2)-kf(:,i+1)*H]'*fai'*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'*R)*e(:,i+1);
end 
for i=2:Bushu
    mxjian(:,i+1)=xjian(:,i)+pp(:,n*(i-2)+1:n*(i-1))*[eye(2)-kf(:,i)*H]'*fai'*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'*R)*e(:,i+1);
    mmxjian(:,i+1)=mxjian(:,i)+pp(:,n*(i-2)+1:n*(i-1))*[eye(2)-kf(:,i)*H]'*fai'*[eye(2)-kf(:,i+2)*H]'*fai'*H'*inv(H*pp(:,n*(i-1)+1:n*i)*H'*R)*e(:,i+1);
end

%-------------------��ͼ����--------------------%
t=1:Bushu;
subplot(2,2,1);plot(t,x(1,t),'b',t,xjian(1,t),'r:');legend('�˲���','��ֵ');%axis([10 100 -150 150]);legend('����Ԥ����');
subplot(2,2,2);plot(t,x(2,t),'b',t,xjian(2,t),'r:');legend('�˲���','��ֵ');%axis([10 100 -150 150]);
subplot(2,2,3);plot(t,x(1,t),'b',t,xxxjian(1,t),'r:');legend('����Ԥ����','��ֵ');
subplot(2,2,4);plot(t,x(2,t),'b',t,xxxjian(2,t),'r:');legend('����Ԥ����','��ֵ');
figure
subplot(2,2,1);plot(t,x(1,t),'b',t,mmxjian(1,t),'r:');legend('����ƽ����','��ֵ');
subplot(2,2,2);plot(t,x(2,t),'b',t,mmxjian(2,t),'r:');legend('����ƽ����','��ֵ');

