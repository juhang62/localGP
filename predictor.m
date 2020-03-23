function [f, v] = predictor(x_star)

global ModelInfo
x_u = ModelInfo.x_u;
u = ModelInfo.u;
y= u;
mu=meanfun(x_star); %prior mean; ad hoc, need more flexible

hyp = ModelInfo.hyp;
%D = size(x_u,2);
D=1; %ad hoc

K1 = knn(x_star, x_u, hyp(1:end-1),0); %same x_star takes a row

psi = K1;

L=ModelInfo.L;
Kinvy=(L'\(L\(y-meanfun(x_u))));
f = mu + psi*Kinvy;

alpha = (L'\(L\psi')); %note the tranpose of psi

v = knn(x_star, x_star, ModelInfo.hyp(1:end-1),0) - psi*alpha;

