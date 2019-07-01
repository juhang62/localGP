function out=kkp_normconst(xh,V,yh,W,m,S,hyp) %can handle all input at once, ie., xh 3*101
VW=V+W;
VWinv=inv(VW);
ntr=size(xh,2);

% tic
% A=VWinv*V;
% B=VWinv*W;
% out=zeros(ntr);
% for r=1:ntr %can vectorize??
%     for s=1:ntr
%         q=VWinv*(V*xh(:,r)+W*yh(:,s));
%         qmm=q-m;
%         term1=exp(hyp(1))^2*exp(-0.5*(xh(:,r)'*V*xh(:,r)+yh(:,s)'*W*yh(:,s)-q'*VW*q)) ...
%             /sqrt(det(S*(VW)+eye(size(S))));
%         term2=exp(-0.5*qmm'*((S+VWinv)\qmm));
%         out(r,s)=term1*term2;
%     end
% end
% toc

SiVWinv=S-S*(VWinv+S)*S;
SiVWinvVW=SiVWinv*VW;
Sm=SiVWinvVW*m;
Vx=V*xh; Wy=W*yh;
test=repmat(dot(xh, Vx)',1,ntr)+repmat(dot(yh, Wy),ntr,1)...
    -repmat(dot(Vx,SiVWinv*Vx)',1,ntr)-repmat(dot(Wy, SiVWinv*Wy),ntr,1) ...
    -2*Vx'*SiVWinv*Wy...
    +2*repmat((Vx'*(Sm-m)),1,ntr)...
    +2*repmat((Wy'*(Sm-m))',ntr,1)...
    +m'*(VW-VW*SiVWinvVW)*m;
haha=exp(hyp(1))^2/sqrt(det(S*(VW)+eye(size(S))))*exp(-0.5*test);
% if any(any((abs(haha-out)>1e-10)))
%     keyboard
% end
out=haha;

end