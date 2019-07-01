function K=DkDx(x, y, hyp, i)

logsigma = hyp(1);
logtheta = hyp(2);

n_x = size(x,1);
n_y = size(y,1);

x = x*ones(1,n_y);
y = ones(n_x,1)*y';

if i==0
    
    K = 2.*exp(1).^logtheta.*pi.^(-1).*(1+exp(1).^logsigma+exp(1) ...
        .^logtheta.*y.^2).*((1+exp(1).^logsigma+exp(1).^logtheta.*x.^2).*( ...
        1+exp(1).^logsigma+exp(1).^logtheta.*y.^2)).^(-3/2).*(y+exp(1) ...
        .^logsigma.*((-1).*x+y)).*(1+(-1).*(1+exp(1).^logsigma+exp(1) ...
        .^logtheta.*x.^2).^(-1).*(exp(1).^logsigma+exp(1).^logtheta.*x.*y) ...
        .^2.*(1+exp(1).^logsigma+exp(1).^logtheta.*y.^2).^(-1)).^(-1/2);
    
elseif i== 1
    
    K = exp(1).^(logsigma+logtheta).*pi.^(-1).*(1+exp(1).^logsigma+exp(1) ...
        .^logtheta.*x.^2).^(-2).*((1+exp(1).^logsigma+exp(1).^logtheta.* ...
        x.^2).*(1+exp(1).^logsigma+exp(1).^logtheta.*y.^2)).^(-1/2).*(1+ ...
        2.*exp(1).^logsigma+exp(1).^logtheta.*(x.^2+exp(1).^logsigma.*(x+( ...
        -1).*y).^2+y.^2)).^(-1).*(1+(-1).*(1+exp(1).^logsigma+exp(1) ...
        .^logtheta.*x.^2).^(-1).*(exp(1).^logsigma+exp(1).^logtheta.*x.*y) ...
        .^2.*(1+exp(1).^logsigma+exp(1).^logtheta.*y.^2).^(-1)).^(-1/2).*( ...
        exp(1).^(2.*logsigma+logtheta).*(x+(-1).*y).^3+(-1).*exp(1) ...
        .^logtheta.*(4.*x.^3+x.^2.*y+y.^3+exp(1).^logtheta.*x.^2.*(x+(-1) ...
        .*y).*(2.*x.^2+x.*y+y.^2))+exp(1).^logsigma.*((-1).*(2+exp(1) ...
        .^logtheta.*(x+(-1).*y).^2).*(x+exp(1).^logtheta.*x.^2.*(x+(-1).* ...
        y)+2.*y)+(-4).*y.*cosh(logsigma)+4.*x.*sinh(logsigma)));
    
elseif i== 2
    
    K = (-1).*exp(1).^logtheta.*pi.^(-1).*(1+exp(1).^logsigma+exp(1) ...
        .^logtheta.*x.^2).^(-2).*(exp(1).^logsigma.*(x+(-1).*y)+(-1).*y).* ...
        ((1+exp(1).^logsigma+exp(1).^logtheta.*x.^2).*(1+exp(1).^logsigma+ ...
        exp(1).^logtheta.*y.^2)).^(-1/2).*(1+2.*exp(1).^logsigma+exp(1) ...
        .^logtheta.*(x.^2+exp(1).^logsigma.*(x+(-1).*y).^2+y.^2)).^(-1).*( ...
        2+2.*exp(1).^logsigma.*(3+2.*exp(1).^logsigma)+exp(1).^logtheta.*( ...
        1+exp(1).^logsigma).*(x.^2+exp(1).^logsigma.*(x+(-1).*y).^2+y.^2)+ ...
        (-1).*exp(1).^(2.*logtheta).*x.^2.*(x.^2+exp(1).^logsigma.*(x+(-1) ...
        .*y).^2+y.^2)).*(1+(-1).*(1+exp(1).^logsigma+exp(1).^logtheta.* ...
        x.^2).^(-1).*(exp(1).^logsigma+exp(1).^logtheta.*x.*y).^2.*(1+exp( ...
        1).^logsigma+exp(1).^logtheta.*y.^2).^(-1)).^(-1/2);
    
end

end