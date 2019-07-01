function [NLML,D_NLML]=likelihood(hyp)

global ModelInfo
x_u = ModelInfo.x_u;
u = ModelInfo.u;
y=u;
mu=meanfun(x_u); %prior mean; ad hoc, need more flexible
%jitter = ModelInfo.jitter;

sigma_n = exp(hyp(end));

[n, ~] = size(x_u);

Knn = knn(x_u, x_u, hyp(1:end-1), 0);

K = Knn;

K = K + eye(n).*sigma_n;

% Cholesky factorisation
[L,p]=chol(K,'lower');

ModelInfo.L = L;

if p > 0
    fprintf(1,'Covariance is ill-conditioned\n');
end

alpha = L'\(L\(y-mu)); %(K)^(-1)*y
NLML = 0.5*(y-mu)'*alpha + sum(log(diag(L))) + log(2*pi)*n/2;

D_NLML = 0*hyp;
Q =  L'\(L\eye(n)) - alpha*alpha';
for i=1:length(hyp)-1
    DKnn = knn(x_u, x_u, hyp(1:end-1),i);     
    DK = DKnn;
    D_NLML(i) = sum(sum(Q.*DK))/2; %sym identity 
end

D_NLML(end) = sigma_n*trace(Q)/2;


