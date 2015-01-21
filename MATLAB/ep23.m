clc;
clear all;
N=9000;
randn('seed',14);
w=randn(6,N);
y(1:6,1)=w(:,1);
for i=2:N
y(:,i)=y(:,i-1)+w(:,i);
end
subplot(2,1,1)
t=1:N;
plot(t,y(:,t))
