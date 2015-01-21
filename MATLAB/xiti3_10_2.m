clc;clear;Bushu=200;
% * * * * * * ����ϵͳģ�� * * * * * * %
Q=1;R=1;fai=0.5;gama=1;H=1;
randn('seed',4);w=sqrt(Q)*randn(1,Bushu+10);
randn('seed',1);v=sqrt(R)*randn(1,Bushu+10);
x(:,1)=0;
y(1)=H*x(:,1)+v(1);
for i=2:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
% * * * * * * ����Riccati���� * * * * * * %
n=1; %ά��
pp(:,:)=10^5*eye(n);
for i=2:Bushu+1;
    temp=pp(:,n*(i-2)+1:n*(i-1));
    kf(:,i-1)=temp*H'*inv(H*temp*H'+R);
    pp(:,n*(i-1)+1:n*i)=fai*[temp-temp*H'*inv(H*temp*H'+R)*H*temp]*fai'+gama*Q*gama';
    kp(:,i-1)=fai*kf(:,i-1);
end
 % * * * * * * ��̬�˲���Ԥ�� * * * * * * %
 xjian(:,1)=zeros(n,1);xxjian(:,1)=zeros(n,1);e(1)=y(1)-H*xxjian(:,1);
 pusaif=(eye(n)-kf(:,Bushu)*H)*fai;
 pusaip=fai*(eye(n)-kf(:,Bushu)*H);
 Qe=H*temp*H'+R;
 for i=2:Bushu+10
     xjian(:,i)=pusaif*xjian(:,i-1)+kf(:,Bushu)*y(i);
     e(i-1)=y(i-1)-H*xxjian(:,i-1);
     xxjian(:,i)=fai*xxjian(:,i-1)+kp(:,Bushu)*e(i-1);
     xxnjian(:,i)=fai*xxjian(:,i); %2��Ԥ��
 end
 % * * * * * * �����������ֵ�� * * * * * * %
for i=1:Bushu+5
    M(1,i)=Q*gama'*H'*inv(Qe);
    M(2,i)=Q*gama'*pusaip'*H'*inv(Qe);
end
for i=1:Bushu
    wjian(1,i)=M(1,i+1)*e(i+1);             %һ��ƽ��
    wjian(2,i)=wjian(1,i)+M(2,i+2)*e(i+2);  %����ƽ��
end
% * * * * * * �۲��������ֵ�� * * * * * * %
for i=1:Bushu+5
    Mv(1,i)=R*inv(Qe);
    Mv(2,i)=-R*kf(:,Bushu)'*fai'*H'*inv(Qe);
    Mv(3,i)=-R*kf(:,Bushu)'*fai'*pusaip'*H'*inv(Qe);
end
for i=1:Bushu
    vjian(1,i)=Mv(1,i)*e(i);                 %һ��ƽ��
    vjian(2,i)=vjian(1,i)+Mv(2,i+1)*e(i+1);  %����ƽ��
    vjian(3,i)=vjian(2,i)+Mv(3,i+2)*e(i+2);  %����ƽ��
end
 figure
 t=1:Bushu;
 subplot(2,2,1);plot(t,x(1,t),t,xjian(1,t),'r:');  %��̬Kalman�˲���
 subplot(2,2,2);plot(t,x(1,t),t,xxnjian(1,t),'r:');%��̬Kalman**2��**Ԥ����
 subplot(2,2,3);plot(t,w(1,t),t,wjian(1,t),'r:');  %���������1��Ԥ��
 subplot(2,2,4);plot(t,w(1,t),t,wjian(2,t),'r:');  %���������2��Ԥ��
 figure
 t=1:Bushu;
 subplot(1,2,1);plot(t,v(1,t),t,vjian(1,t),'r:');  %�۲������0��Ԥ��
 subplot(1,2,2);plot(t,v(1,t),t,vjian(2,t),'r:');  %�۲������1��Ԥ��
  figure
 subplot(1,1,1);plot(t,v(1,t),t,vjian(3,t),'r:');  %�۲������2��Ԥ��
 