%¸½Â¼I ÎÈÌ¬KalmanÂË²¨Ëã·¨Maltlab·ÂÕæÍ¨Ê½
function[xigema,kf,kp,pusaif,pusaip,e,xjian,xxjian]=TS_WTKalman(fai,gama,H,Q,R,y)
%-----------------------------------%

%-----------------------------------%
n=size(fai,1);
L=size(R,1);
Bushu=size(y,2);
pp(:,:)=eye(n);
%-----------------------------------%
for i=2:Bushu+1;
    temp=pp(:,n*(i-2)+1:n*(i-1));
    pp(:,n*(i-1)+1:n*i)=fai*[temp-temp*H'*inv(H*temp*H'+R)*H*temp]*fai'+gama*Q*gama';
    kf(:,L*(i-1)+1:L*(i-1)+L)=temp*H'*inv(H*temp*H'+R);
    kp(:,L*(i-1)+1:L*(i-1)+L)=fai*kf(:,L*(i-1)+1:L*(i-1)+L);
end
 %---------------------------------------%
 kf=kf(:,end-L+1:end);kp=kp(:,end-L+1:end);
 pusaif=(eye(n)-kf*H)*fai;pusaip=fai*(eye(n)-kf*H);
 xigema=pp(:,end-n+1:end);
 
 xjian(:,Bushu)=zeros(1,n)';
 for i=2:Bushu
     xjian(:,i)=pusaif*xjian(:,i-1)+kf*y(:,i);
 end
 %---------------------------------------%

 xxjian(:,Bushu)=zeros(1,n)';
 for i=2:Bushu
     e(:,i)=y(:,i)-H*xxjian(:,i-1);
     xxjian(:,i)=pusaip*xxjian(:,i-1)+kp*y(:,i);
 end

