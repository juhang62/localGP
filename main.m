
%close all;
clear
addpath ./Utilities
addpath ./Kernels

rng(1)

global ModelInfo

set(0,'defaulttextinterpreter','latex')

%% Setup
jitter = 1e-8;
ModelInfo.jitter=jitter;
% sig_noi = 0.01;

%load('stochastic fisher2.mat')
load('bistable2d.mat')
ModelInfo.spdim=2;
ModelInfo.sparse=1;
ModelInfo.ord=0;
ModelInfo.nx=length(xgrid);
ndata=size(Utr,1);
t0=20; %better be 50 for 2d
xtrain=formx12d(Utr(:,1),ModelInfo.spdim,ModelInfo.ord);
%ytrain=U(:,t0+1)+normrnd(0,sig_noi,ndata,1);
ytrain=Utr(:,2);

ModelInfo.x_u = xtrain;
ModelInfo.u = ytrain;


ModelInfo.hyp = log([0.01 .01 .01 exp(-6)]); %logged parameter: front sig, 2 length scales in denominator,obs noise; squared built in


[ModelInfo.hyp,~,~] = minimize(ModelInfo.hyp, @likelihood, -25);
%ModelInfo.hyp=[-8.0183/2   -1.4782   -1.4782/2   -9.2253]; %1d 
%ModelInfo.hyp=[ -3.7380    2.0935    1.5359   -9.1735]; %2d
[NLML,~]=likelihood(ModelInfo.hyp);

%xstar = formx2d(U(:,t0+1)); %formated input
%[Kpred, Kvar] = predictor(xstar);
%S=eye(length(xstar))*sig_noi^2;
% if ModelInfo.sparse==1
%     S=spalloc(ndata,ndata,10*ndata);
%     S=S+speye(ndata)*sig_noi^2;
% else
%     S=eye(ndata)*sig_noi^2;
% end
totcell=25;
ucell=cell(totcell,1);
Scell=cell(totcell,1);
% [unew, Snew] = predictorGau(U(:,t0+1),S); %make sure Snew is p.s.d
[unew, Snew] = predictor(formx12d(Utr(:,2),ModelInfo.spdim,ModelInfo.ord)); %make sure Snew is p.s.d
Snew=Snew+exp(ModelInfo.hyp(end))*eye(size(Snew));
ucell{1}=unew;
Scell{1}=Snew;
%plotstuff(xgrid,unew,diag(Snew))
for i=1:totcell-1
    i
[unew, Snew] = predictorGau(unew,Snew);
ucell{i+1}=unew;
Scell{i+1}=Snew;
%figure
%plotstuff(xgrid,unew,diag(Snew))
%keyboard
end
Umean=mean(Ute,3);

