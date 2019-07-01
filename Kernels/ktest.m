function out=ktest(x, y, hyp, i)

logsigma = hyp(1); %front sigma
logtheta1 = hyp(2); %length scale
logtheta2 = hyp(3); %length scale
% n_x = size(x,1);
% n_y = size(y,1);
% 
% x = x*ones(1,n_y);
% y = ones(n_x,1)*y';
[ndata,dim_x]=size(x);
difsq1=zeros(ndata,ndata);
for j=1:dim_x-1
    x1=x(:,j);
    y1=y(:,j);
    difsq1=difsq1+(x1-y1').^2;
end

    x1=x(:,end);
    y1=y(:,end);
    difsq2=(x1-y1').^2;



if i==0 || i==1
    
    
    %%%squared Exp 
    out = exp(1).^(logsigma+(-1/2).*exp(1).^((-1).*logtheta1).*difsq1+(-1/2).*exp(1).^((-1).*logtheta2).*difsq2);
    %out = exp(1).^(logsigma+(-1/2).*exp(1).^((-1).*logtheta1).*difsq);

% 
%     %%%%squared exp
%     out = exp(1).^(logsigma+(-1/2).*exp(1).^((-1).*logtheta).*difsq); %gradient wrt logsigma
    
    
elseif i== 2
    %keyboard
%     K = exp(1).^logtheta.*pi.^(-1).*((1+exp(1).^logsigma+exp(1) ...
%         .^logtheta.*x.^2).*(1+exp(1).^logsigma+exp(1).^logtheta.*y.^2)).^( ...
%         -3/2).*(1+(-1).*(1+exp(1).^logsigma+exp(1).^logtheta.*x.^2).^(-1) ...
%         .*(exp(1).^logsigma+exp(1).^logtheta.*x.*y).^2.*(1+exp(1) ...
%         .^logsigma+exp(1).^logtheta.*y.^2).^(-1)).^(-1/2).*((-1).*exp(1) ...
%         .^(2.*logsigma).*(x+(-1).*y).^2+exp(1).^(logsigma+logtheta).*x.*( ...
%         x+(-1).*y).^2.*y+(-1).*exp(1).^logsigma.*(x.^2+(-4).*x.*y+y.^2)+ ...
%         x.*y.*(2+exp(1).^logtheta.*(x.^2+y.^2)));
%     
    %%%squared exp
    %out = (1/2).*exp(1).^(logsigma+(-1).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*difsq).*difsq;
    
        out = (1/2).*exp(1).^(logsigma+(-1).*logtheta1+(-1/2).*exp(1).^((-1).* ...
        logtheta1).*difsq1+(-1/2).*exp(1).^((-1).*logtheta2).*difsq2).*difsq1;
    
    %out = (1/2).*exp(1).^(logsigma+(-1).*logtheta1+(-1/2).*exp(1).^((-1).*logtheta).*difsq).*difsq;
    
elseif i==3
        out = (1/2).*exp(1).^(logsigma+(-1).*logtheta2+(-1/2).*exp(1).^((-1).* ...
        logtheta1).*difsq1+(-1/2).*exp(1).^((-1).*logtheta2).*difsq2).*difsq2;
%out=zeros(ndata,ndata);

end

end