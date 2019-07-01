function out=k(x, y, hyp, i)

logsigma = hyp(1); %front sigma
logtheta = hyp(2); %length scale

% n_x = size(x,1);
% n_y = size(y,1);
% 
% x = x*ones(1,n_y);
% y = ones(n_x,1)*y';
[ndata,dim_x]=size(x);
difsq=zeros(ndata,ndata);
for j=1:dim_x
    x1=x(:,j);
    y1=y(:,j);
    difsq=difsq+(x1-y1').^2;
end


if i==0
    
    %K = 2.*pi.^(-1).*asin((exp(1).^logsigma+exp(1).^logtheta.*x.*y).*((1+exp(1).^logsigma+exp(1).^logtheta.*x.^2).*(1+exp(1).^logsigma+exp(1).^logtheta.*y.^2)).^(-1/2));
    
    %%%squared Exp 
    out = exp(1).^(logsigma+(-1/2).*exp(1).^((-1).*logtheta).*difsq);
    
elseif i== 1
    %keyboard
%     K = exp(1).^logsigma.*pi.^(-1).*((1+exp(1).^logsigma+exp(1) ...
%         .^logtheta.*x.^2).*(1+exp(1).^logsigma+exp(1).^logtheta.*y.^2)).^( ...
%         -3/2).*(2+2.*exp(1).^logsigma+exp(1).^(logsigma+logtheta).*(x+(-1) ...
%         .*y).^2+(-1).*exp(1).^(2.*logtheta).*x.*(x+(-1).*y).^2.*y+2.*exp( ...
%         1).^logtheta.*(x.^2+(-1).*x.*y+y.^2)).*(1+(-1).*(1+exp(1) ...
%         .^logsigma+exp(1).^logtheta.*x.^2).^(-1).*(exp(1).^logsigma+exp(1) ...
%         .^logtheta.*x.*y).^2.*(1+exp(1).^logsigma+exp(1).^logtheta.*y.^2) ...
%         .^(-1)).^(-1/2);

    %%%%squared exp
    out = exp(1).^(logsigma+(-1/2).*exp(1).^((-1).*logtheta).*difsq); %gradient wrt logsigma
    
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
    out = (1/2).*exp(1).^(logsigma+(-1).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*difsq).*difsq;

end

end