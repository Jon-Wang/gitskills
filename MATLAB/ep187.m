%例3.8.2 多道ARMA信号wiener滤波器仿真源程序
clc;clear all;Bushu=300;
%---------------生成系统模型----------------%
A1=[-1 0;0.9 0.2];C1=[0.95 -1.1;-1.4 0.6];Q=eye(2);R=4*eye(2);
randn('seed',1);w=sqrt(Q)*randn(2,Bushu+10);
randn('seed',2);v=sqrt(R)*randn(2,Bushu+10);
s(:,1)=[0;0];y(:,1)=s(:,1)+v(:,1);
for i=2:Bushu+10
    s(:,i)=-A1*s(:,i-1)+C1*w(:,i-1);y(:,i)=s(:,i)+v(:,i);
end
fai=-A1;gama=C1;H=eye(2);
%-------------Riccati方程---------------%
n=2;
pp(:,:)=eye(n);
 for i=2:Bushu+1;
    temp=pp(:,n*(i-2)+1:n*(i-1));
    pp(:,n*(i-1)+1:n*i)=fai*[temp-temp*H'*inv(H*temp*H'+R)*H*temp]*fai'+gama*Q*gama';
    kf(:,:,i)=temp*H'*inv(H*temp*H'+R);
 end
 kf=kf(:,end-1:end);
 kp=fai*kf;pusaip=fai*(eye(n)-kf*H);xigema=pp(:,end-1:end);
 %-------------解耦部分---------------%
 syms q;
 digits(4);
 f=eye(2)-pusaip*q;
 phisym=det(f);
 phipoly=sym2poly(phisym);
 adj=simple((inv(f)*det(f)));
 
 Lambda=simple(vpa(phisym*eye(2)-H*adj*kp*q));
 Lambda11=sym2poly(Lambda(1,1)); Lambda12=sym2poly(Lambda(1,2));
 Lambda21=sym2poly(Lambda(2,1)); Lambda22=sym2poly(Lambda(2,2));
 Lambda0=[Lambda11(3) Lambda12(3);Lambda21(3) Lambda22(3)];
 Lambda1=[Lambda11(2) Lambda12(2);Lambda21(2) Lambda22(2)];
 Lambda2=[Lambda11(1) Lambda12(1);Lambda21(1) Lambda22(1)];
 
 J=simple(vpa(phisym*eye(2)*q));
 J11=sym2poly(J(1,1));J12=sym2poly(J(1,2));
 J21=sym2poly(J(2,1));J22=sym2poly(J(2,2));
 
 J0=[J11(4) 0;0 J22(4)]; J1=[J11(3) 0;0 J22(3)];
 J2=[J11(2) 0;0 J22(2)]; J3=[J11(1) 0;0 J22(1)];
 
 Mv0=R*inv(H*xigema*H'+R); Mv1=-R*kf'*fai'*H'*inv(H*xigema*H'+R);
 Ln0=Mv1; Ln1=Mv0;
 
 LnLam0=Mv1*Lambda0;LnLam1=Mv1*Lambda1+Mv0*Lambda0;
 LnLam2=Mv1*Lambda2+Mv0*Lambda1;LnLam3=Mv0*Lambda2;
 
 K10=J0-LnLam0; K11=J1-LnLam1; K12=J2-LnLam2; K13=J3-LnLam3;
 sjian(:,1:2)=zeros(2);
 for i=3:Bushu
     sjian(:,i)=-phipoly(2)*sjian(:,i-1)-phipoly(1)*sjian(:,i-2)+K10*y(:,i+1)+K11*y(:,i)+K12*y(:,i-1)+K13*y(:,i-2);
 end
 
 %-------------作图部分------------%
 t=1:Bushu;
 subplot(2,2,1);plot(t,s(1,t),'b',t,sjian(1,t),'r:');
 subplot(2,2,2);plot(t,s(2,t),'b',t,sjian(2,t),'r:');

