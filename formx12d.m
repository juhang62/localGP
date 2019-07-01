function out=formx12d(u,spdim,ord)
[n1,n2]=size(u);
if n2==1
    nx=sqrt(n1);
elseif n1==n2 
    nx=n1;
elseif n1~=n2
    keyboard
end


if spdim==2
    ran_m=reshape(1:nx*nx,[nx nx]);
    ran1 = [ran_m(:,end) ran_m(:,1:end-1)];
    ran2 = [ran_m(end,:) ; ran_m(1:end-1,:)];
    ran3 = [ran_m(2:end,:); ran_m(1,:)];
    ran4 = [ran_m(:,2:end) ran_m(:,1)];
    out=[u(ran1(:)) u(ran2(:)) u(ran3(:)) u(ran4(:)) u(:) ];
else
    col1=[u(end);u(1:end-1)];
    col3=[u(2:end);u(1)];
    out=[col1 col3 u];
end
    
    
if ord==1    
    out(:,1:end-1)=sort(out(:,1:end-1),2);
end

%ran1(abs(1-neib(i))+1)=nx*nx;
end

% 
% %%%%%%%%%%%neib 2d unord
%   2
% 1 5 4
%   3
