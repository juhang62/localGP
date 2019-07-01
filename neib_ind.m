function indout=neib_ind(ind,nx, spdim, ord, m)
if spdim==2
[s1,s2]=ind2sub([nx nx],ind);
neibsub=repmat([s1 s2], [5,1]);
ad=[0 -1; -1 0; 1 0; 0 1; 0 0];
neibsub=neibsub+ad;
neibsub(neibsub==0)=nx;
neibsub(neibsub==nx+1)=1;
indout=sub2ind([nx nx],neibsub(:,1),neibsub(:,2));
else
    if ind==1
        indout=[nx,ind+1,ind]';
    elseif ind==nx
        indout=[ind-1,1,ind]';
    else
        indout=[ind-1,ind+1,ind]';
    end
end

if ord==1
    mj=m(indout);
    [~,sortind]=sort(mj(1:end-1));
    tmp=indout(1:end-1);
    indout(1:end-1)=tmp(sortind);
end
end