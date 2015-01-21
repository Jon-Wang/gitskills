%��2.9.1 ����RLS�㷨�������
clear;clc;Bushu=1000;
Q=[1 0;0 0.81];A=[0.9 0;0.7 0.5];
randn('seed',1);e=sqrt(Q)*randn(2,Bushu);
y(:,1)=e(:,1);
for i=2:Bushu
    y(:,i)=A*y(:,i-1)+e(:,i);%����AR()ģ��
end
%-----------------ģ�Ͳ�����ʶ------------------%
Na=1;%ģ�͵Ĵ�������Na
N=2;%Ԥ��������ά��N
Num=N*Na;
p0=10^5*eye(Num);
fai(:,1)=[0,0]';p(:,1:2)=10^5*eye(2);
sita(1:4,1)=zeros(4,1);%��Ϊ����A����4��δ֪����
for k=1:N
    for i=2:Bushu
        fai(:,i)=[y(1,i-1),y(2,i-1)]';
        sita((k-1)*Num+1:k*Num,i)=sita((k-1)*Num+1:k*Num,i-1)+p(:,Num*(i-2)+1:Num*(i-1))*fai(:,i)/(1+fai(:,i)'*p(:,Num*(i-2)+1:Num*(i-1))*fai(:,i))*(y(k,i)-fai(:,i)'*sita((k-1)*Num+1:k*Num,i-1));
        p(:,Num*(i-1)+1:Num*i)=p(:,Num*(i-2)+1:Num*(i-1))-p(:,Num*(i-2)+1:Num*(i-1))*fai(:,i)*fai(:,i)'*p(:,Num*(i-2)+1:Num*(i-1))/(1+fai(:,i)'*p(:,Num*(i-2)+1:Num*(i-1))*fai(:,i));
    end
end  
%----------------����������Ĺ���----------------%
for k=1:N
    for i=1:Bushu
        emixiu(i)=y(k,i)-fai(:,i)'*sita((k-1)*Num+1:k*Num,i);
    end
    taoe(k,1)=emixiu(1)^2;
    for i=2:Bushu
        taoe(k,i)=taoe(k,i-1)+1/i*[emixiu(i)^2-taoe(k,i-1)];
    end
end    
%--------------------��ͼ����--------------------%
t=1:Bushu;
subplot(2,2,1);plot(t,sita(1,t),'r',t,sita(2,t),'r',t,sita(3,t),'r',t,sita(4,t),'r');
axis([0 Bushu -0.5 1.5]);
line([0,Bushu],[0.9,0.9]);line([0,Bushu],[0,0]);
line([0,Bushu],[0.7,0.7]);line([0,Bushu],[0.5,0.5]);
subplot(2,2,2);plot(t,taoe(1,t),'r',t,taoe(2,t),'r');
line([0,Bushu],[1,1]);line([0,Bushu],[0.81,0.81]);

