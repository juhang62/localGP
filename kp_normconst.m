function [const, varargout]=kp_normconst(xh,V,m,S,hyp) %now handle multiple training data xh: 3*ndata
invS=inv(S);
term1=exp(hyp(1))/sqrt(det(S*V+eye(size(V))));

%ntr=size(xh,2);
%const=zeros(1,ntr);
% for r=1:ntr
% xm=xh(:,r)-m;
% Vm=V*xm;
% term2=exp(-0.5*(xm'*(Vm-V*((invS+V)\Vm))));
% const(r)=term1*term2;
% end

xm=xh-m;
Vm=V*xm;
term2=exp(-0.5*(dot(xm,(Vm-V*((invS+V)\Vm)))));
const=term1*term2;

if nargout==2
    tmp_mean=((V+invS)\(V*xh+invS*m));
    varargout{1}= const.*tmp_mean(end,:);
end
end