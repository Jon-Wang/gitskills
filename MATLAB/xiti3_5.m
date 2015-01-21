clear all;clc;Bushu=300;
fai=0.8;gama=1;H=1;R=0.1;
n=1;%*****1维
%%%%%%%*********Bernoulli-Gaussian白噪声**********%%%%%
Taog=25;Lmd=0.15;Q=Lmd*Taog;
randn('seed',13);g=sqrt(Taog)*randn(1,Bushu+10);
rand('state',1);k=rand(1,Bushu+10);
for i=1:Bushu+10
    if k(i)<=Lmd
        b(i)=1;
    else
        b(i)=0;
    end
    w(i)=b(i)*g(i);
end
%%%%%%%%**********生成V(t)**********
randn('seed',3);v=sqrt(R)*randn(1,Bushu+10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(:,1)=0;
y(1,1)=H*x(:,1)+v(1);
for i=2:Bushu+10
    x(:,i)=fai*x(:,i-1)+gama*w(i-1);
    y(i)=H*x(:,i)+v(i);
end
%%%%%%%%%%*********Kalman滤波器估值部分***********
e(1)=y(1);k(:,1)=0;p(:,:)=0;xjian(:,1)=0;
for i=1:Bushu+8
    xxjian(:,i+1)=fai*xjian(:,i);
    e(:,i+1)=y(i+1)-H*xxjian(:,i+1);
    pp(:,i)=fai*p(:,i)*fai'+gama*Q*gama';
    kf(:,i+1)=pp(:,i)*H'*inv(H*pp(:,i)*H'+R);
    xjian(:,i+1)=xxjian(:,i+1)+kf(:,i+1)*e(:,i+1);
    p(:,(i+1))=[1-kf(:,i+1)*H]*pp(:,i);
end
%%%%%%%%%%%*********白噪声估值器部分**********
N=3;%%%%%%%%3步平滑
for i=1:Bushu+8
    pusaip(:,:,i)=fai*(eye(n)-kf(:,i)*H);
    Qe(:,:,i)=H*pp(:,i)*H'+R;
end
for i=1:Bushu+5
    M(1,i)=Q*gama'*H'*inv(Qe(:,:,i+1));
    M(2,i)=Q*gama'*pusaip(:,:,i+1)'*H'*inv(Qe(:,:,i+2));
    M(3,i)=Q*gama'*pusaip(:,:,i+2)'*pusaip(:,:,i+1)'*H'*inv(Qe(:,:,i+3));
end
for i=1:Bushu
    wjian(1,i)=M(1,i+1)*e(i+1);                  %%%%%一步平滑
    wjian(2,i)=wjian(1,i)+M(2,i+2)*e(i+2);       %%%%%二步平滑
    wjian(3,i)=wjian(2,i)+M(3,i+3)*e(i+3);       %%%%%三步平滑
end
%%%%%%%*********plot*******
for tt=1:N
    subplot(2,2,tt)
    t=1:Bushu;plot(t,wjian(tt,t),'r.')
    for t=1:Bushu
        hh=line([t,t],[0,w(t)]);
        set(hh,'color','k')
    end
end

    


