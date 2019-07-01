function [f, v] = predictorGau(m,S)

global ModelInfo
x_u = ModelInfo.x_u;
u = ModelInfo.u;
y= u;
[ndata,ndim]=size(x_u);
hyp = ModelInfo.hyp;
L=ModelInfo.L;
Kinvy=(L'\(L\(y-meanfun(x_u)))); %only about training data; beta in GP note

%Lnn=lnn(x_u,m_xform,S,hyp);

%[Eukb, v, BB, CC]=deal(zeros(ndata));
%[Eukb, v]=deal(zeros(ndata));
if ModelInfo.sparse==1
  [Eukb, v]=deal(spalloc(ndata,ndata,10*ndata));
else
  [Eukb, v]=deal(zeros(ndata));
end
Lam=eye(ndim)*exp(hyp(2));
invLam=eye(ndim)*(1/exp(hyp(2)));

if length(hyp)==4
    Lam(end,end)=exp(hyp(3));
    invLam(end,end)=1/exp(hyp(3));
end
%%%%%%%%%%
Lnn=zeros(ndata);
for j=1:ndata
%     j
%     pause(0.1)
    nind=neib_ind(j,ModelInfo.nx,ModelInfo.spdim,ModelInfo.ord,m);
    mj=m(nind);
    Sb=S(nind, nind);
    Lnn(:,j)=kp_normconst(x_u',invLam,mj,Sb,hyp);
end

%%%%%%%%%%%
parfor j=1:ndata
    %keyboard
    nind=neib_ind(j,ModelInfo.nx,ModelInfo.spdim,ModelInfo.ord,m);
    mj=m(nind);
    Sb=S(nind, nind);
    Lj=chol(Lam+Sb,'lower');
    tmp_mean=(mj-Sb*(Lj'\(Lj\mj)))+(x_u'-Lam*(Lj'\(Lj\x_u')));   %ndim*ndata
    %     Eukb_slice=zeros(5,1);
    %     v_slice=zeros(5,1);
    %tic
    for i=1:ndata
        %for i=nind'
        ind= nind==i;
        if ModelInfo.sparse==1 && ~any(ind)
            continue
        end
        
        if any(ind) %can combine to kp?
            Eukb(i,j)=(Lnn(:,j).*(tmp_mean(ind,:))')'*Kinvy;
        else
            nindp1=[nind;i];
            Sbi=S(nindp1,nindp1);
            mji=[mj; m(i)];
            PEuk=[eye(ndim) zeros(ndim,1)];
            [~,Euk]=kp_normconst(PEuk'*x_u',PEuk'*invLam*PEuk,mji,Sbi,hyp);
            Eukb(i,j)=Euk*Kinvy; %need normalize!
            %keyboard
        end
        %%%%%%%%%%%
        %[i j]
        if i<=j %|| 1==1  %%!!!!!!!!!!!!!!!!!
            nindi=neib_ind(i,ModelInfo.nx,ModelInfo.spdim,ModelInfo.ord,m);
            [nindij,ia,ic]=unique([nindi; nind],'stable');
            P=[eye(ndim) zeros(ndim,length(nindij)-ndim)];
            Q=zeros(size(P));
            subtmp=ic(ndim+1:end);
            indtmp=sub2ind(size(Q),(1:ndim)',subtmp);
            Q(indtmp)=1;
            P_Ek=P-Q;
            mij=m(nindij);
            Sbij=S(nindij,nindij);
            %                 if any(P*nindij~=nindi)
            %                     keyboard
            %                 end
            
            
            %[i j]
            Lrs=kkp_normconst(P'*x_u', P'*invLam*P, Q'*x_u', Q'*invLam*Q, mij, Sbij, hyp);
            bstuffb=Kinvy'*(Lrs-Lnn(:,i)*(Lnn(:,j))')*Kinvy; %not sym! need fix
            ETr=trace(L'\(L\Lrs));
            Ek=kp_normconst(0,P_Ek'*invLam*P_Ek,mij,Sbij,hyp);
            v(i,j)=S(i,j)+bstuffb-(m(i)*Lnn(:,j)'+m(j)*Lnn(:,i)')*Kinvy+Ek-ETr;
            %             BB(i,j)=Ek;
            %             CC(i,j)=ETr;
%             if i==j
%                 keyboard
%             end
        end
    end
    %toc
end
%%%%%%%%%%%%%%%%%%

v=v+v'-diag(diag(v));
v=v+Eukb+Eukb';
f = m + Lnn'*Kinvy;
if any(diag(v)<0)
error('fuck')
end
% keyboard
end


